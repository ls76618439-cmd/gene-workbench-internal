from __future__ import annotations

import copy
import difflib
import hashlib
import os
import sys
import traceback
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from mcp.server.fastmcp import FastMCP
from pydna.amplify import pcr as pydna_pcr
from pydna.assembly2 import (
    gibson_assembly as pydna_gibson_assembly,
    golden_gate_assembly as pydna_golden_gate_assembly,
    ligation_assembly as pydna_ligation_assembly,
    restriction_ligation_assembly as pydna_restriction_ligation_assembly,
)
from pydna.dseqrecord import Dseqrecord
from pydna.design import assembly_fragments as pydna_assembly_fragments, primer_design as pydna_primer_design
from pydna.primer import Primer
from primer3.bindings import calc_hairpin, calc_homodimer, calc_heterodimer, calc_tm, design_primers as primer3_design_primers

APP_HOME = Path(os.environ.get("LOCALAPPDATA", str(Path.home()))) / "GeneWorkbench"
DATA = APP_HOME / "data"
OUTPUTS = APP_HOME / "outputs"
DATA.mkdir(parents=True, exist_ok=True)
OUTPUTS.mkdir(parents=True, exist_ok=True)

mcp = FastMCP("gene-workbench")
VERSION = "1.2.0"
IUPAC_DNA = set("ACGTRYSWKMBDHVN")


def _format_for(path: Path) -> str:
    ext = path.suffix.lower()
    if ext in {".gb", ".gbk", ".genbank"}:
        return "genbank"
    if ext in {".fa", ".fasta", ".fna", ".fas"}:
        return "fasta"
    if ext == ".dna":
        return "snapgene"
    raise ValueError(f"Unsupported sequence file: {path.name}")


def _read_file(path: Path) -> SeqRecord:
    record = SeqIO.read(str(path), _format_for(path))
    if "molecule_type" not in record.annotations:
        record.annotations["molecule_type"] = "DNA"
    return record


def _record_path(sequence_id: str) -> Path:
    path = DATA / f"{sequence_id}.gb"
    if not path.exists():
        raise FileNotFoundError(f"Unknown sequence_id: {sequence_id}")
    return path


def _load(sequence_id: str) -> SeqRecord:
    return SeqIO.read(str(_record_path(sequence_id)), "genbank")


def _clean_dna(text: str, allow_empty: bool = False) -> str:
    value = text.upper().replace(" ", "").replace("\n", "").replace("\r", "").replace("\t", "")
    if not value and not allow_empty:
        raise ValueError("DNA sequence is empty")
    if any(base not in IUPAC_DNA for base in value):
        raise ValueError("sequence must contain IUPAC DNA characters only")
    return value


def _save_record(record: SeqRecord, label: str) -> tuple[str, Path]:
    seq_text = str(record.seq).upper()
    digest = hashlib.sha256(seq_text.encode()).hexdigest()[:16]
    safe_label = "".join(c if c.isalnum() or c in "-_" else "_" for c in label).strip("_") or "sequence"
    sequence_id = f"{safe_label}_{digest}"
    path = DATA / f"{sequence_id}.gb"
    record.id = sequence_id
    record.name = sequence_id[:16]
    record.annotations["molecule_type"] = "DNA"
    SeqIO.write(record, str(path), "genbank")
    return sequence_id, path


def _new_record_like(record: SeqRecord, seq: Seq) -> SeqRecord:
    out = SeqRecord(seq, id=record.id, name=record.name, description=record.description)
    out.annotations = copy.deepcopy(record.annotations)
    out.annotations["molecule_type"] = "DNA"
    out.dbxrefs = list(record.dbxrefs)
    return out


def _remap_simple_location(
    loc: SimpleLocation, start0: int, end0: int, replacement_len: int
) -> SimpleLocation | None:
    s = int(loc.start)
    e = int(loc.end)
    delta = replacement_len - (end0 - start0)

    if e <= start0:
        ns, ne = s, e
    elif s >= end0:
        ns, ne = s + delta, e + delta
    elif s < start0 and e > end0:
        ns, ne = s, e + delta
    elif s < start0 < e <= end0:
        ns, ne = s, start0
    elif start0 <= s < end0 < e:
        ns, ne = start0 + replacement_len, e + delta
    else:
        return None

    if ne <= ns:
        return None
    return SimpleLocation(ns, ne, strand=loc.strand)


def _remap_location_edit(
    location: Any, start0: int, end0: int, replacement_len: int
) -> Any | None:
    if isinstance(location, CompoundLocation):
        parts = []
        for part in location.parts:
            mapped = _remap_simple_location(part, start0, end0, replacement_len)
            if mapped is not None:
                parts.append(mapped)
        if not parts:
            return None
        if len(parts) == 1:
            return parts[0]
        return CompoundLocation(parts, operator=location.operator)
    return _remap_simple_location(location, start0, end0, replacement_len)


def _feature_intersects_edit(location: Any, start0: int, end0: int) -> bool:
    parts = location.parts if isinstance(location, CompoundLocation) else [location]
    if start0 == end0:
        return any(int(part.start) < start0 < int(part.end) for part in parts)
    return any(int(part.start) < end0 and int(part.end) > start0 for part in parts)


def _edit_record(record: SeqRecord, start0: int, end0: int, replacement: str) -> tuple[SeqRecord, dict[str, int]]:
    replacement_seq = Seq(replacement)
    new_seq = record.seq[:start0] + replacement_seq + record.seq[end0:]
    out = _new_record_like(record, new_seq)
    remapped = 0
    dropped = 0
    overlapping = 0
    translation_qualifiers_removed = 0
    for feat in record.features:
        intersects = _feature_intersects_edit(feat.location, start0, end0)
        mapped_location = _remap_location_edit(feat.location, start0, end0, len(replacement))
        if mapped_location is None:
            dropped += 1
            continue
        new_feat = copy.deepcopy(feat)
        new_feat.location = mapped_location
        if intersects:
            overlapping += 1
            if "translation" in new_feat.qualifiers:
                new_feat.qualifiers.pop("translation", None)
                translation_qualifiers_removed += 1
            notes = list(new_feat.qualifiers.get("note", []))
            notes.append("GeneWorkbench: feature overlaps an edited interval; sequence-derived qualifiers require revalidation.")
            new_feat.qualifiers["note"] = notes
        out.features.append(new_feat)
        remapped += 1
    return out, {
        "features_remapped": remapped,
        "features_dropped": dropped,
        "features_overlapping_edit": overlapping,
        "translation_qualifiers_removed": translation_qualifiers_removed,
    }


def _to_dseqrecord(record: SeqRecord) -> Dseqrecord:
    circular = str(record.annotations.get("topology", "")).lower() == "circular"
    return Dseqrecord.from_SeqRecord(record, circular=circular)


def _record_from_dseq(product: Dseqrecord, description: str) -> SeqRecord:
    seq_text = str(product.seq.watson)
    out = SeqRecord(Seq(seq_text), id="product", name="product", description=description)
    out.annotations = copy.deepcopy(getattr(product, "annotations", {}))
    out.annotations["molecule_type"] = "DNA"
    out.annotations["topology"] = "circular" if bool(getattr(product, "circular", False)) else "linear"
    out.features = copy.deepcopy(getattr(product, "features", []))
    return out


def _save_assembly_products(products: list[Dseqrecord], label: str, max_products: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for i, product in enumerate(products[:max_products], start=1):
        record = _record_from_dseq(product, f"{label} assembly product {i}")
        sequence_id, path = _save_record(record, f"{label}_{i}")
        if sequence_id in seen:
            continue
        seen.add(sequence_id)
        rows.append(
            {
                "sequence_id": sequence_id,
                "length_bp": len(record.seq),
                "topology": record.annotations.get("topology"),
                "feature_count": len(record.features),
                "stored_path": str(path),
            }
        )
    return rows


def _enzyme_objects(enzymes: list[str]) -> list[Any]:
    if not enzymes:
        raise ValueError("at least one restriction enzyme is required")
    batch = RestrictionBatch(enzymes)
    return list(batch)


def _is_circular(record: SeqRecord) -> bool:
    return str(record.annotations.get("topology", "")).lower() == "circular"


def _enzyme_profile(enzyme: Any) -> dict[str, Any]:
    fst5 = getattr(enzyme, "fst5", None)
    fst3 = getattr(enzyme, "fst3", None)
    scd5 = getattr(enzyme, "scd5", None)
    scd3 = getattr(enzyme, "scd3", None)
    size = int(getattr(enzyme, "size", len(str(getattr(enzyme, "site", "")))))
    one_cut_pair = scd5 is None and scd3 is None
    primary_cut_outside = fst5 is not None and (int(fst5) < 0 or int(fst5) > size)
    type_iis_geometry = bool(one_cut_pair and primary_cut_outside)
    return {
        "name": str(enzyme),
        "recognition_site": str(enzyme.site),
        "recognition_length": size,
        "fst5": fst5,
        "fst3": fst3,
        "overhang_length": abs(int(getattr(enzyme, "ovhg", 0))),
        "blunt": bool(enzyme.is_blunt()),
        "cuts_outside_recognition_site": primary_cut_outside,
        "single_cut_pair": one_cut_pair,
        "golden_gate_type_iis_geometry": type_iis_geometry,
        "elucidate": enzyme.elucidate(),
    }


def _golden_gate_enzyme_objects(enzymes: list[str], allow_blunt: bool) -> tuple[list[Any], list[dict[str, Any]]]:
    objs = _enzyme_objects(enzymes)
    profiles = [_enzyme_profile(e) for e in objs]
    incompatible = [
        p["name"]
        for p in profiles
        if not p["golden_gate_type_iis_geometry"] or (p["blunt"] and not allow_blunt)
    ]
    if incompatible:
        raise ValueError(
            "assembly_golden_gate requires enzymes with a single cut pair outside the recognition site"
            + (" and non-blunt ends" if not allow_blunt else "")
            + f"; incompatible: {', '.join(incompatible)}. Use assembly_ligation for ordinary restriction/ligation."
        )
    return objs, profiles


def _primer_metrics(full_sequence: str, annealing_sequence: str) -> dict[str, Any]:
    full_sequence = str(full_sequence).upper()
    annealing_sequence = str(annealing_sequence).upper()
    hp = calc_hairpin(full_sequence)
    hd = calc_homodimer(full_sequence)
    return {
        "full_sequence_5to3": full_sequence,
        "annealing_sequence_5to3": annealing_sequence,
        "tail_sequence_5to3": full_sequence[: max(0, len(full_sequence) - len(annealing_sequence))],
        "full_length": len(full_sequence),
        "annealing_length": len(annealing_sequence),
        "annealing_tm_c": round(float(calc_tm(annealing_sequence)), 2),
        "hairpin_tm_c": round(float(hp.tm), 2),
        "homodimer_tm_c": round(float(hd.tm), 2),
    }


def _mapped_location_from_boundaries(location: Any, boundary_map: dict[int, int]) -> Any | None:
    def one(loc: SimpleLocation) -> SimpleLocation | None:
        s = int(loc.start)
        e = int(loc.end)
        if s not in boundary_map or e not in boundary_map:
            return None
        ns = boundary_map[s]
        ne = boundary_map[e]
        if ne <= ns:
            return None
        return SimpleLocation(ns, ne, strand=loc.strand)

    if isinstance(location, CompoundLocation):
        parts = [one(p) for p in location.parts]
        if any(p is None for p in parts):
            return None
        clean = [p for p in parts if p is not None]
        if len(clean) == 1:
            return clean[0]
        return CompoundLocation(clean, operator=location.operator)
    return one(location)


@mcp.tool()
def sequence_import(path: str, label: str = "sequence") -> dict[str, Any]:
    """Import a local FASTA, GenBank, or SnapGene .dna file into the local sequence store."""
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(path)
    record = _read_file(source)
    sequence_id, stored = _save_record(record, label)
    return {
        "sequence_id": sequence_id,
        "length_bp": len(record.seq),
        "topology": record.annotations.get("topology"),
        "features": len(record.features),
        "stored_path": str(stored),
        "source_name": source.name,
    }


@mcp.tool()
def sequence_info(sequence_id: str) -> dict[str, Any]:
    """Return sequence metadata without returning the full sequence."""
    record = _load(sequence_id)
    seq = str(record.seq).upper()
    gc = 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
    return {
        "sequence_id": sequence_id,
        "length_bp": len(seq),
        "gc_percent": round(gc, 3),
        "topology": record.annotations.get("topology"),
        "description": record.description,
        "feature_count": len(record.features),
        "sha256": hashlib.sha256(seq.encode()).hexdigest(),
    }


@mcp.tool()
def sequence_extract(sequence_id: str, start: int, end: int, strand: int = 1) -> dict[str, Any]:
    """Extract a 1-based inclusive interval. strand may be 1 or -1."""
    if start < 1 or end < start:
        raise ValueError("Use 1-based inclusive coordinates with end >= start")
    record = _load(sequence_id)
    if end > len(record.seq):
        raise ValueError("end exceeds sequence length")
    subseq = record.seq[start - 1 : end]
    if strand == -1:
        subseq = subseq.reverse_complement()
    elif strand != 1:
        raise ValueError("strand must be 1 or -1")
    return {"sequence_id": sequence_id, "start": start, "end": end, "strand": strand, "sequence": str(subseq)}


@mcp.tool()
def sequence_find(sequence_id: str, query: str, both_strands: bool = True, max_hits: int = 100) -> dict[str, Any]:
    """Find exact nucleotide occurrences and return 1-based coordinates."""
    record = _load(sequence_id)
    target = str(record.seq).upper()
    q = _clean_dna(query)
    hits: list[dict[str, int]] = []
    pos = target.find(q)
    while pos >= 0 and len(hits) < max_hits:
        hits.append({"start": pos + 1, "end": pos + len(q), "strand": 1})
        pos = target.find(q, pos + 1)
    if both_strands:
        rc = str(Seq(q).reverse_complement())
        if rc != q:
            pos = target.find(rc)
            while pos >= 0 and len(hits) < max_hits:
                hits.append({"start": pos + 1, "end": pos + len(rc), "strand": -1})
                pos = target.find(rc, pos + 1)
    hits.sort(key=lambda x: (x["start"], -x["strand"]))
    return {"sequence_id": sequence_id, "query_length": len(q), "hit_count": len(hits), "hits": hits}


@mcp.tool()
def sequence_features(sequence_id: str, feature_type: str | None = None, max_features: int = 500) -> dict[str, Any]:
    """List annotated features with coordinates and selected qualifiers."""
    record = _load(sequence_id)
    rows = []
    for feat in record.features:
        if feature_type and feat.type.lower() != feature_type.lower():
            continue
        rows.append(
            {
                "type": feat.type,
                "start": int(feat.location.start) + 1,
                "end": int(feat.location.end),
                "strand": feat.location.strand,
                "label": (feat.qualifiers.get("label") or feat.qualifiers.get("gene") or feat.qualifiers.get("locus_tag") or [None])[0],
                "note": (feat.qualifiers.get("note") or [None])[0],
            }
        )
        if len(rows) >= max_features:
            break
    return {"sequence_id": sequence_id, "count": len(rows), "features": rows}


@mcp.tool()
def sequence_replace(sequence_id: str, start: int, end: int, replacement: str, label: str = "edited") -> dict[str, Any]:
    """Replace a 1-based inclusive interval, remap unaffected annotations, and save a new GenBank record."""
    record = _load(sequence_id)
    if start < 1 or end < start or end > len(record.seq):
        raise ValueError("invalid coordinates")
    replacement = _clean_dna(replacement)
    new_record, stats = _edit_record(record, start - 1, end, replacement)
    sequence_id_new, path = _save_record(new_record, label)
    return {
        "source_sequence_id": sequence_id,
        "sequence_id": sequence_id_new,
        "old_length_bp": len(record.seq),
        "new_length_bp": len(new_record.seq),
        "replaced_start": start,
        "replaced_end": end,
        "replacement_length": len(replacement),
        **stats,
        "stored_path": str(path),
    }


@mcp.tool()
def sequence_insert(sequence_id: str, position: int, insert_sequence: str, label: str = "inserted") -> dict[str, Any]:
    """Insert DNA before a 1-based position. Use position=len+1 to append. Features are remapped."""
    record = _load(sequence_id)
    if position < 1 or position > len(record.seq) + 1:
        raise ValueError("position must be between 1 and sequence length + 1")
    insert_sequence = _clean_dna(insert_sequence)
    start0 = position - 1
    new_record, stats = _edit_record(record, start0, start0, insert_sequence)
    sequence_id_new, path = _save_record(new_record, label)
    return {
        "source_sequence_id": sequence_id,
        "sequence_id": sequence_id_new,
        "position": position,
        "insert_length": len(insert_sequence),
        "old_length_bp": len(record.seq),
        "new_length_bp": len(new_record.seq),
        **stats,
        "stored_path": str(path),
    }


@mcp.tool()
def sequence_delete(sequence_id: str, start: int, end: int, label: str = "deleted") -> dict[str, Any]:
    """Delete a 1-based inclusive interval, remap unaffected annotations, and save a new record."""
    record = _load(sequence_id)
    if start < 1 or end < start or end > len(record.seq):
        raise ValueError("invalid coordinates")
    new_record, stats = _edit_record(record, start - 1, end, "")
    sequence_id_new, path = _save_record(new_record, label)
    return {
        "source_sequence_id": sequence_id,
        "sequence_id": sequence_id_new,
        "deleted_start": start,
        "deleted_end": end,
        "deleted_length": end - start + 1,
        "old_length_bp": len(record.seq),
        "new_length_bp": len(new_record.seq),
        **stats,
        "stored_path": str(path),
    }


@mcp.tool()
def sequence_diff(source_sequence_id: str, target_sequence_id: str, max_changes: int = 100) -> dict[str, Any]:
    """Compare two stored sequences in their stored coordinate systems and summarize substitutions, insertions, and deletions."""
    source_record = _load(source_sequence_id)
    target_record = _load(target_sequence_id)
    source = str(source_record.seq).upper()
    target = str(target_record.seq).upper()
    matcher = difflib.SequenceMatcher(None, source, target, autojunk=False)
    changes = []
    counts = {"replace": 0, "insert": 0, "delete": 0}
    truncated = False
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue
        if len(changes) >= max_changes:
            truncated = True
            break
        counts[tag] = counts.get(tag, 0) + 1
        old_fragment = source[i1:i2]
        new_fragment = target[j1:j2]
        changes.append(
            {
                "type": tag,
                "source_start": i1 + 1,
                "source_end": i2,
                "target_start": j1 + 1,
                "target_end": j2,
                "source_length": len(old_fragment),
                "target_length": len(new_fragment),
                "source_sequence": old_fragment if len(old_fragment) <= 200 else old_fragment[:200] + "...",
                "target_sequence": new_fragment if len(new_fragment) <= 200 else new_fragment[:200] + "...",
            }
        )
    return {
        "source_sequence_id": source_sequence_id,
        "target_sequence_id": target_sequence_id,
        "source_length_bp": len(source),
        "target_length_bp": len(target),
        "change_count": len(changes),
        "change_type_counts": counts,
        "truncated": truncated,
        "changes": changes,
        "source_topology": source_record.annotations.get("topology"),
        "target_topology": target_record.annotations.get("topology"),
        "scope": "Coordinate-based sequence comparison. Circular sequences are not origin-normalized, so a pure origin rotation can appear as a large diff.",
    }


@mcp.tool()
def feature_remap(source_sequence_id: str, target_sequence_id: str, label: str = "remapped") -> dict[str, Any]:
    """Transfer annotations onto a related target sequence using unchanged coordinate blocks; mark features whose sequence changed."""
    source_record = _load(source_sequence_id)
    target_record = _load(target_sequence_id)
    source = str(source_record.seq).upper()
    target = str(target_record.seq).upper()
    matcher = difflib.SequenceMatcher(None, source, target, autojunk=False)
    boundary_map: dict[int, int] = {}
    for a, b, size in matcher.get_matching_blocks():
        for offset in range(size + 1):
            boundary_map[a + offset] = b + offset

    out = _new_record_like(target_record, target_record.seq)
    mapped = 0
    dropped = 0
    sequence_changed = 0
    translation_qualifiers_removed = 0
    for feat in source_record.features:
        new_location = _mapped_location_from_boundaries(feat.location, boundary_map)
        if new_location is None:
            dropped += 1
            continue
        new_feat = copy.deepcopy(feat)
        source_feature_seq = str(feat.extract(source_record.seq)).upper()
        new_feat.location = new_location
        target_feature_seq = str(new_feat.extract(target_record.seq)).upper()
        if source_feature_seq != target_feature_seq:
            sequence_changed += 1
            if "translation" in new_feat.qualifiers:
                new_feat.qualifiers.pop("translation", None)
                translation_qualifiers_removed += 1
            notes = list(new_feat.qualifiers.get("note", []))
            notes.append("GeneWorkbench: mapped feature sequence differs from source; sequence-derived qualifiers require revalidation.")
            new_feat.qualifiers["note"] = notes
        out.features.append(new_feat)
        mapped += 1

    sequence_id_new, path = _save_record(out, label)
    return {
        "source_sequence_id": source_sequence_id,
        "target_sequence_id": target_sequence_id,
        "sequence_id": sequence_id_new,
        "features_mapped": mapped,
        "features_dropped": dropped,
        "features_sequence_changed": sequence_changed,
        "translation_qualifiers_removed": translation_qualifiers_removed,
        "stored_path": str(path),
        "scope": "Best-effort annotation transfer using unchanged sequence blocks. Repetitive or rearranged sequences can produce ambiguous mappings and should be reviewed.",
    }


@mcp.tool()
def sequence_translate(sequence_id: str, start: int, end: int, strand: int = 1, table: int = 1, to_stop: bool = False) -> dict[str, Any]:
    """Translate a nucleotide interval using an NCBI translation table."""
    piece = sequence_extract(sequence_id, start, end, strand)["sequence"]
    protein = str(Seq(piece).translate(table=table, to_stop=to_stop))
    return {"sequence_id": sequence_id, "start": start, "end": end, "strand": strand, "protein": protein, "aa_length": len(protein)}


@mcp.tool()
def restriction_analyze(sequence_id: str, enzymes: list[str]) -> dict[str, Any]:
    """Find restriction enzyme cut sites, respecting linear versus circular topology."""
    record = _load(sequence_id)
    batch = RestrictionBatch(enzymes)
    linear = not _is_circular(record)
    analysis = batch.search(record.seq, linear=linear)
    return {
        "sequence_id": sequence_id,
        "topology": "circular" if not linear else "linear",
        "enzymes": {str(e): [int(x) for x in pos] for e, pos in analysis.items()},
        "enzyme_profiles": [_enzyme_profile(e) for e in batch],
    }


@mcp.tool()
def primer_design(sequence_id: str, target_start: int, target_length: int, product_min: int = 120, product_max: int = 1200, num_return: int = 5) -> dict[str, Any]:
    """Design candidate primer pairs around a 1-based target using Primer3 and return Primer3 quality metrics."""
    record = _load(sequence_id)
    template = str(record.seq).upper()
    if target_start < 1 or target_length < 1 or target_start + target_length - 1 > len(template):
        raise ValueError("target must be a valid 1-based interval inside the template")
    if product_min < 1 or product_max < product_min:
        raise ValueError("product size range is invalid")
    args = {"SEQUENCE_ID": sequence_id, "SEQUENCE_TEMPLATE": template, "SEQUENCE_TARGET": [target_start - 1, target_length]}
    global_args = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 57.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 35.0,
        "PRIMER_MAX_GC": 65.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
        "PRIMER_NUM_RETURN": num_return,
    }
    result = primer3_design_primers(args, global_args)
    pairs = []
    for i in range(int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))):
        pairs.append(
            {
                "left": result.get(f"PRIMER_LEFT_{i}_SEQUENCE"),
                "right": result.get(f"PRIMER_RIGHT_{i}_SEQUENCE"),
                "primer3_left_coordinate": result.get(f"PRIMER_LEFT_{i}"),
                "primer3_right_coordinate": result.get(f"PRIMER_RIGHT_{i}"),
                "left_tm": result.get(f"PRIMER_LEFT_{i}_TM"),
                "right_tm": result.get(f"PRIMER_RIGHT_{i}_TM"),
                "left_gc_percent": result.get(f"PRIMER_LEFT_{i}_GC_PERCENT"),
                "right_gc_percent": result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT"),
                "left_self_any_th": result.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
                "right_self_any_th": result.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
                "left_self_end_th": result.get(f"PRIMER_LEFT_{i}_SELF_END_TH"),
                "right_self_end_th": result.get(f"PRIMER_RIGHT_{i}_SELF_END_TH"),
                "left_hairpin_th": result.get(f"PRIMER_LEFT_{i}_HAIRPIN_TH"),
                "right_hairpin_th": result.get(f"PRIMER_RIGHT_{i}_HAIRPIN_TH"),
                "pair_compl_any_th": result.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
                "pair_compl_end_th": result.get(f"PRIMER_PAIR_{i}_COMPL_END_TH"),
                "pair_penalty": result.get(f"PRIMER_PAIR_{i}_PENALTY"),
                "product_size": result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
            }
        )
    return {
        "sequence_id": sequence_id,
        "pair_count": len(pairs),
        "pairs": pairs,
        "primer3_explain": {
            "left": result.get("PRIMER_LEFT_EXPLAIN"),
            "right": result.get("PRIMER_RIGHT_EXPLAIN"),
            "pair": result.get("PRIMER_PAIR_EXPLAIN"),
        },
        "scope": "Primer3 candidate design on the provided template; off-target specificity outside this template is not assessed.",
    }


@mcp.tool()
def pcr_simulate(
    sequence_id: str,
    forward_primer: str,
    reverse_primer: str,
    anneal_min: int = 13,
    label: str = "pcr_product",
) -> dict[str, Any]:
    """Simulate PCR with 5'-to-3' forward and reverse primers and save the product."""
    if anneal_min < 8:
        raise ValueError("anneal_min must be at least 8")
    record = _load(sequence_id)
    fwd = _clean_dna(forward_primer)
    rev = _clean_dna(reverse_primer)
    product = pydna_pcr(Primer(fwd), Primer(rev), _to_dseqrecord(record), limit=anneal_min)
    out = _record_from_dseq(product, f"PCR product from {sequence_id}")
    sequence_id_new, path = _save_record(out, label)
    return {
        "template_sequence_id": sequence_id,
        "sequence_id": sequence_id_new,
        "product_length_bp": len(out.seq),
        "topology": out.annotations.get("topology"),
        "stored_path": str(path),
        "scope": "In-silico PCR on the supplied template and primer sequences. This does not establish wet-lab amplification efficiency, specificity outside the supplied template, or absence of secondary-structure effects.",
    }


@mcp.tool()
def assembly_gibson(
    sequence_ids: list[str],
    overlap_min: int = 25,
    circular_only: bool = True,
    max_products: int = 10,
    label: str = "gibson",
) -> dict[str, Any]:
    """Assemble fragments with Gibson-style terminal homology and save candidate products."""
    if len(sequence_ids) < 2:
        raise ValueError("Gibson assembly requires at least two fragments")
    if overlap_min < 10:
        raise ValueError("overlap_min must be at least 10")
    try:
        fragments = [_to_dseqrecord(_load(sid)) for sid in sequence_ids]
        products = pydna_gibson_assembly(fragments, limit=overlap_min, circular_only=circular_only)
        rows = _save_assembly_products(products, label, max_products)
    except BaseException as exc:
        raise RuntimeError("Gibson assembly failed:\n" + traceback.format_exc()) from exc
    return {
        "input_sequence_ids": sequence_ids,
        "product_count": len(rows),
        "circular_only": circular_only,
        "overlap_min": overlap_min,
        "products": rows,
    }


@mcp.tool()
def golden_gate_assess(sequence_ids: list[str], enzymes: list[str]) -> dict[str, Any]:
    """Report enzyme cut geometry and recognition-site counts for supplied sequences without selecting a preferred enzyme."""
    if not sequence_ids:
        raise ValueError("at least one sequence is required")
    objs = _enzyme_objects(enzymes)
    profiles = [_enzyme_profile(e) for e in objs]
    sequence_rows = []
    for sid in sequence_ids:
        record = _load(sid)
        linear = not _is_circular(record)
        analysis = RestrictionBatch(enzymes).search(record.seq, linear=linear)
        sequence_rows.append(
            {
                "sequence_id": sid,
                "topology": "circular" if not linear else "linear",
                "sites": {str(e): [int(x) for x in pos] for e, pos in analysis.items()},
                "site_counts": {str(e): len(pos) for e, pos in analysis.items()},
            }
        )
    return {
        "sequence_ids": sequence_ids,
        "enzyme_profiles": profiles,
        "sequences": sequence_rows,
        "scope": "Factual recognition/cut geometry and site counts only; this tool does not choose an enzyme or claim experimental success.",
    }


@mcp.tool()
def gibson_primer_design(
    sequence_ids: list[str],
    overlap: int = 30,
    target_tm: float = 60.0,
    circular: bool = True,
    maxlink: int = 40,
) -> dict[str, Any]:
    """Design candidate PCR primers with 5' homology tails for Gibson-style assembly in the supplied fragment order."""
    import warnings

    if len(sequence_ids) < 2:
        raise ValueError("at least two fragments are required")
    if overlap < 15:
        raise ValueError("overlap must be at least 15 bp")
    if maxlink < 0:
        raise ValueError("maxlink must be non-negative")

    source_records = [_load(sid) for sid in sequence_ids]
    input_topologies = [
        "circular" if _is_circular(record) else "linear"
        for record in source_records
    ]
    warning_messages: list[str] = []
    if any(topology == "circular" for topology in input_topologies):
        warning_messages.append(
            "One or more supplied templates are circular. Primer design treats each sequence_id as the exact fragment/template to amplify end-to-end; verify that these are the intended fragment boundaries."
        )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        base_amplicons = [pydna_primer_design(_to_dseqrecord(record), target_tm=target_tm) for record in source_records]
        designed = pydna_assembly_fragments(base_amplicons, overlap=overlap, maxlink=maxlink, circular=circular)
        warning_messages.extend(str(w.message) for w in caught)

    rows = []
    for sid, source_record, amp in zip(sequence_ids, source_records, designed):
        f_full = str(amp.forward_primer.seq).upper()
        r_full = str(amp.reverse_primer.seq).upper()
        f_anneal = str(amp.forward_primer.footprint).upper()
        r_anneal = str(amp.reverse_primer.footprint).upper()
        simulated = pydna_pcr(Primer(f_full), Primer(r_full), _to_dseqrecord(source_record), limit=13)
        pair_dimer = calc_heterodimer(f_full, r_full)
        rows.append(
            {
                "source_sequence_id": sid,
                "forward_primer": _primer_metrics(f_full, f_anneal),
                "reverse_primer": _primer_metrics(r_full, r_anneal),
                "primer_pair_heterodimer_tm_c": round(float(pair_dimer.tm), 2),
                "simulated_pcr_product_length_bp": len(simulated),
            }
        )
    return {
        "input_sequence_ids": sequence_ids,
        "fragment_order": sequence_ids,
        "input_topologies": input_topologies,
        "input_semantics": "Each sequence_id is treated as the exact fragment/template intended for end-to-end PCR amplification before assembly.",
        "circular": circular,
        "requested_overlap_bp": overlap,
        "target_annealing_tm_c": target_tm,
        "designs": rows,
        "warnings": warning_messages,
        "scope": "Candidate overlap-tailed primers computed for the supplied fragment order. PCR simulation uses only the supplied template; genome-wide or host-wide off-target specificity and wet-lab efficiency are not assessed.",
    }


@mcp.tool()
def assembly_golden_gate(
    sequence_ids: list[str],
    enzymes: list[str],
    allow_blunt: bool = False,
    circular_only: bool = True,
    max_products: int = 10,
    label: str = "golden_gate",
) -> dict[str, Any]:
    """Run Golden Gate-style restriction/ligation using enzymes that cut outside their recognition site."""
    if not sequence_ids:
        raise ValueError("at least one fragment is required")
    enzyme_objects, profiles = _golden_gate_enzyme_objects(enzymes, allow_blunt)
    fragments = [_to_dseqrecord(_load(sid)) for sid in sequence_ids]
    products = pydna_golden_gate_assembly(
        fragments,
        enzyme_objects,
        allow_blunt=allow_blunt,
        circular_only=circular_only,
    )
    rows = _save_assembly_products(products, label, max_products)
    return {
        "input_sequence_ids": sequence_ids,
        "enzymes": enzymes,
        "enzyme_profiles": profiles,
        "product_count": len(rows),
        "circular_only": circular_only,
        "products": rows,
        "scope": "Computational complete-digest/ligation assembly from the supplied sequences and enzyme definitions.",
    }


@mcp.tool()
def assembly_ligation(
    sequence_ids: list[str],
    enzymes: list[str] | None = None,
    allow_blunt: bool = True,
    circular_only: bool = True,
    max_products: int = 10,
    label: str = "ligation",
) -> dict[str, Any]:
    """Assemble DNA by ligation. If enzymes are given, digest-and-ligate is simulated first."""
    if not sequence_ids:
        raise ValueError("at least one fragment is required")
    fragments = [_to_dseqrecord(_load(sid)) for sid in sequence_ids]
    if enzymes:
        products = pydna_restriction_ligation_assembly(
            fragments,
            _enzyme_objects(enzymes),
            allow_blunt=allow_blunt,
            circular_only=circular_only,
        )
        mode = "restriction_ligation"
    else:
        products = pydna_ligation_assembly(
            fragments,
            allow_blunt=allow_blunt,
            circular_only=circular_only,
        )
        mode = "ligation"
    rows = _save_assembly_products(products, label, max_products)
    return {
        "input_sequence_ids": sequence_ids,
        "mode": mode,
        "enzymes": enzymes or [],
        "product_count": len(rows),
        "circular_only": circular_only,
        "products": rows,
    }


@mcp.tool()
def construct_validate(sequence_id: str, max_issues: int = 100) -> dict[str, Any]:
    """Run computational sequence, annotation, and CDS consistency checks without claiming experimental validation."""
    record = _load(sequence_id)
    seq = str(record.seq).upper()
    issues: list[dict[str, Any]] = []

    invalid_chars = sorted(set(seq) - IUPAC_DNA)
    if invalid_chars:
        issues.append({"severity": "error", "type": "invalid_sequence_characters", "characters": invalid_chars})

    feature_bounds_ok = True
    for index, feat in enumerate(record.features):
        s0 = int(feat.location.start)
        e0 = int(feat.location.end)
        if s0 < 0 or e0 > len(record.seq) or e0 <= s0:
            feature_bounds_ok = False
            issues.append(
                {
                    "severity": "error",
                    "type": "feature_bounds",
                    "feature_index": index,
                    "feature_type": feat.type,
                    "start": s0 + 1,
                    "end": e0,
                }
            )
            if len(issues) >= max_issues:
                break

    cds_checks = []
    for index, feat in enumerate(record.features):
        if feat.type.lower() != "cds":
            continue
        cds_seq = feat.extract(record.seq)
        codon_start_raw = (feat.qualifiers.get("codon_start") or ["1"])[0]
        try:
            codon_start = int(codon_start_raw)
        except (TypeError, ValueError):
            codon_start = 1
            if len(issues) < max_issues:
                issues.append({"severity": "warning", "type": "invalid_codon_start", "feature_index": index, "value": codon_start_raw})
        if codon_start not in {1, 2, 3}:
            codon_start = 1
            if len(issues) < max_issues:
                issues.append({"severity": "warning", "type": "invalid_codon_start", "feature_index": index, "value": codon_start_raw})

        coding = cds_seq[codon_start - 1 :]
        frame_ok = len(coding) % 3 == 0
        table_raw = (feat.qualifiers.get("transl_table") or ["1"])[0]
        try:
            table = int(table_raw)
            codon_table = CodonTable.unambiguous_dna_by_id[table]
        except (TypeError, ValueError, KeyError):
            table = 1
            codon_table = CodonTable.unambiguous_dna_by_id[1]
            if len(issues) < max_issues:
                issues.append({"severity": "warning", "type": "invalid_translation_table", "feature_index": index, "value": table_raw})

        trimmed = coding[: len(coding) - (len(coding) % 3)] if len(coding) % 3 else coding
        protein = str(trimmed.translate(table=table, to_stop=False)) if trimmed else ""
        internal_stops = protein[:-1].count("*") if protein else 0
        first_codon = str(coding[:3]).upper() if len(coding) >= 3 else None
        last_codon = str(coding[-3:]).upper() if len(coding) >= 3 and frame_ok else None
        start_codon_recognized = first_codon in codon_table.start_codons if first_codon and set(first_codon) <= set("ACGT") else None
        terminal_stop_recognized = last_codon in codon_table.stop_codons if last_codon and set(last_codon) <= set("ACGT") else None

        declared_translation = "".join(feat.qualifiers.get("translation", [])).replace(" ", "").replace("\\n", "").upper() or None
        computed_translation = protein.rstrip("*")
        translation_matches_annotation = None
        interpretation_qualifiers = [
            key
            for key in ("transl_except", "ribosomal_slippage", "pseudo", "pseudogene")
            if feat.qualifiers.get(key) is not None
        ]
        interpretation_limited = bool(interpretation_qualifiers)
        if declared_translation is not None and not interpretation_limited:
            translation_matches_annotation = declared_translation == computed_translation

        label = (feat.qualifiers.get("label") or feat.qualifiers.get("gene") or feat.qualifiers.get("locus_tag") or [None])[0]
        cds_checks.append(
            {
                "feature_index": index,
                "label": label,
                "start": int(feat.location.start) + 1,
                "end": int(feat.location.end),
                "strand": feat.location.strand,
                "translation_table": table,
                "coding_length_bp": len(coding),
                "frame_ok": frame_ok,
                "first_codon": first_codon,
                "start_codon_recognized_for_table": start_codon_recognized,
                "last_codon": last_codon,
                "terminal_stop_recognized_for_table": terminal_stop_recognized,
                "internal_stop_count": internal_stops,
                "declared_translation_present": declared_translation is not None,
                "translation_matches_annotation": translation_matches_annotation,
                "interpretation_limited": interpretation_limited,
                "interpretation_qualifiers": interpretation_qualifiers,
            }
        )
        if not interpretation_limited and not frame_ok and len(issues) < max_issues:
            issues.append({"severity": "warning", "type": "cds_frame", "feature_index": index, "label": label})
        if not interpretation_limited and internal_stops and len(issues) < max_issues:
            issues.append(
                {
                    "severity": "warning",
                    "type": "cds_internal_stop",
                    "feature_index": index,
                    "label": label,
                    "count": internal_stops,
                }
            )
        if translation_matches_annotation is False and len(issues) < max_issues:
            issues.append(
                {
                    "severity": "warning",
                    "type": "cds_translation_annotation_mismatch",
                    "feature_index": index,
                    "label": label,
                }
            )

    errors = [x for x in issues if x.get("severity") == "error"]
    warnings_out = [x for x in issues if x.get("severity") == "warning"]
    return {
        "sequence_id": sequence_id,
        "length_bp": len(record.seq),
        "topology": record.annotations.get("topology"),
        "sequence_characters_valid": not invalid_chars,
        "feature_bounds_valid": feature_bounds_ok,
        "cds_count": len(cds_checks),
        "cds_checks": cds_checks,
        "error_count": len(errors),
        "warning_count": len(warnings_out),
        "structural_checks_pass": len(errors) == 0,
        "computational_checks_pass": len(errors) == 0 and len(warnings_out) == 0,
        "overall_pass": len(errors) == 0 and len(warnings_out) == 0,
        "issues": issues,
        "scope": "Computational consistency checks only. This is not experimental validation of expression, cloning efficiency, biological function, or sequence identity of a physical sample.",
    }


@mcp.tool()
def sequence_export(sequence_id: str, output_name: str, format: str = "genbank") -> dict[str, Any]:
    """Export a stored sequence as GenBank or FASTA into the local GeneWorkbench outputs folder."""
    fmt = format.lower()
    if fmt not in {"genbank", "fasta"}:
        raise ValueError("format must be genbank or fasta")
    record = _load(sequence_id)
    safe_name = Path(output_name).name
    if not safe_name:
        raise ValueError("output_name is empty")
    if fmt == "genbank" and not Path(safe_name).suffix:
        safe_name += ".gb"
    if fmt == "fasta" and not Path(safe_name).suffix:
        safe_name += ".fasta"
    dest = OUTPUTS / safe_name
    SeqIO.write(record, str(dest), fmt)
    return {"sequence_id": sequence_id, "output_path": str(dest), "format": fmt}


@mcp.tool()
def tool_status() -> dict[str, Any]:
    """Report local Gene Workbench status and storage locations."""
    exe_dir = Path(sys.executable).resolve().parent if getattr(sys, "frozen", False) else Path(__file__).resolve().parent
    return {
        "version": VERSION,
        "app_home": str(APP_HOME),
        "data_dir": str(DATA),
        "outputs_dir": str(OUTPUTS),
        "seqkit_installed": (exe_dir / "seqkit.exe").exists(),
        "transport": "stdio",
    }


if __name__ == "__main__":
    mcp.run(transport="stdio")