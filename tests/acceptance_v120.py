import asyncio
import json
import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

ROOT = Path(r"C:\ProgramData\ChatGPT-Operator\workspace\gene_workbench_dist\acceptance_v120")
ROOT.mkdir(parents=True, exist_ok=True)
PYTHON = r"D:\bio-ai-workbench\envs\sequence\Scripts\python.exe"
SERVER = r"D:\gene-workbench-internal\src\gene_workbench.py"


def write_gb(path: Path, seq: str, circular=False, features=None):
    r = SeqRecord(Seq(seq), id=path.stem, name=path.stem, description="acceptance")
    r.annotations["molecule_type"] = "DNA"
    r.annotations["topology"] = "circular" if circular else "linear"
    r.features = features or []
    SeqIO.write(r, str(path), "genbank")


def write_fa(path: Path, seq: str):
    r = SeqRecord(Seq(seq), id=path.stem, description="")
    SeqIO.write(r, str(path), "fasta")


def randseq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def payload(result):
    sc = getattr(result, "structuredContent", None)
    if sc is not None:
        return sc
    sc = getattr(result, "structured_content", None)
    if sc is not None:
        return sc
    for item in getattr(result, "content", []):
        text = getattr(item, "text", None)
        if text:
            try:
                return json.loads(text)
            except Exception:
                return {"text": text}
    return {}


async def call(session, name, args, expect_error=False):
    r = await session.call_tool(name, args)
    is_error = getattr(r, "isError", getattr(r, "is_error", False))
    if expect_error:
        assert is_error, f"{name} unexpectedly succeeded: {r}"
        return r
    if is_error:
        raise RuntimeError(f"{name} failed: {r}")
    return payload(r)


def prepare_inputs():
    # Valid CDS with a matching /translation qualifier; edit crosses the CDS.
    cds = "ATG" + ("GCC" * 98) + "TAA"
    protein = str(Seq(cds).translate()).rstrip("*")
    edit_seq = ("ACGT" * 25) + cds + ("TGCA" * 25)
    features = [
        SeqFeature(
            SimpleLocation(100, 400, strand=1),
            type="CDS",
            qualifiers={"label": ["good_cds"], "translation": [protein], "transl_table": ["1"]},
        ),
        SeqFeature(SimpleLocation(420, 480, strand=1), type="misc_feature", qualifiers={"label": ["tail_feature"]}),
    ]
    write_gb(ROOT / "edit_source.gb", edit_seq, False, features)

    pseudo_seq = "ATGTAGAAAA"
    pseudo_features = [
        SeqFeature(
            SimpleLocation(0, len(pseudo_seq), strand=1),
            type="CDS",
            qualifiers={"label": ["pseudo_cds"], "pseudo": [""]},
        )
    ]
    write_gb(ROOT / "pseudo.gb", pseudo_seq, False, pseudo_features)

    # Circular EcoRI recognition site split across the origin: suffix GAA + prefix TTC.
    write_gb(ROOT / "circular_origin.gb", "TTC" + ("A" * 20) + "GAA", True, [])

    # PCR and Gibson design templates.
    pcr_seq = randseq(600, 7)
    write_gb(ROOT / "pcr_template.gb", pcr_seq, False, [])
    fwd = pcr_seq[:24]
    rev = str(Seq(pcr_seq[-24:]).reverse_complement())

    for i, seed in enumerate((501, 502, 503), start=1):
        write_fa(ROOT / f"design{i}.fa", randseq(500, seed))

    # Gibson assembly fragments with explicit terminal overlaps.
    z = randseq(25, 101)
    x = randseq(25, 102)
    y = randseq(25, 103)
    a = randseq(70, 201)
    b = randseq(70, 202)
    c = randseq(70, 203)
    write_fa(ROOT / "gibson1.fa", z + a + x)
    write_fa(ROOT / "gibson2.fa", x + b + y)
    write_fa(ROOT / "gibson3.fa", y + c + z)

    # Golden Gate fixture adapted from OpenCloning's own assembly test.
    gg = [
        ("gg1.gb", "GGTCTCAattaAAAAAttaaAGAGACC", False),
        ("gg2.gb", "GGTCTCAttaaCCCCCatatAGAGACC", False),
        ("gg3.gb", "GGTCTCAatatGGGGGccggAGAGACC", False),
        ("gg4.gb", "TTTTattaAGAGACCTTTTTGGTCTCAccggTTTT", True),
    ]
    for name, seq, circular in gg:
        write_gb(ROOT / name, seq, circular, [])

    return fwd, rev


async def main():
    fwd, rev = prepare_inputs()
    params = StdioServerParameters(command=PYTHON, args=[SERVER])
    async with stdio_client(params) as (rs, ws):
        async with ClientSession(rs, ws) as s:
            await s.initialize()
            tools = await s.list_tools()
            names = [t.name for t in tools.tools]
            print("TOOL_COUNT", len(names))
            assert len(names) == 22, names

            # Import/edit/annotation safety.
            src = await call(s, "sequence_import", {"path": str(ROOT / "edit_source.gb"), "label": "edit_src"})
            src_id = src["sequence_id"]
            val = await call(s, "construct_validate", {"sequence_id": src_id})
            assert val["computational_checks_pass"] is True
            assert val["cds_checks"][0]["translation_matches_annotation"] is True
            assert val["cds_checks"][0]["start_codon_recognized_for_table"] is True
            assert val["cds_checks"][0]["terminal_stop_recognized_for_table"] is True

            pseudo = await call(s, "sequence_import", {"path": str(ROOT / "pseudo.gb"), "label": "pseudo"})
            pseudo_val = await call(s, "construct_validate", {"sequence_id": pseudo["sequence_id"]})
            assert pseudo_val["cds_checks"][0]["interpretation_limited"] is True
            assert "pseudo" in pseudo_val["cds_checks"][0]["interpretation_qualifiers"]
            pseudo_issue_types = {x["type"] for x in pseudo_val["issues"]}
            assert "cds_frame" not in pseudo_issue_types
            assert "cds_internal_stop" not in pseudo_issue_types

            edited = await call(
                s,
                "sequence_replace",
                {"sequence_id": src_id, "start": 150, "end": 152, "replacement": "AAA", "label": "edited_cds"},
            )
            assert edited["features_overlapping_edit"] >= 1
            assert edited["translation_qualifiers_removed"] == 1
            edited_features = await call(s, "sequence_features", {"sequence_id": edited["sequence_id"], "feature_type": "CDS"})
            assert edited_features["count"] == 1

            # Circular topology must be respected for restriction sites spanning the origin.
            cir = await call(s, "sequence_import", {"path": str(ROOT / "circular_origin.gb"), "label": "circ"})
            digest = await call(s, "restriction_analyze", {"sequence_id": cir["sequence_id"], "enzymes": ["EcoRI"]})
            assert digest["topology"] == "circular"
            assert len(digest["enzymes"]["EcoRI"]) == 1, digest

            # Primer3 must expose quality metrics and scope.
            pcrt = await call(s, "sequence_import", {"path": str(ROOT / "pcr_template.gb"), "label": "pcr_template"})
            primers = await call(
                s,
                "primer_design",
                {"sequence_id": pcrt["sequence_id"], "target_start": 200, "target_length": 60, "product_min": 150, "product_max": 500},
            )
            assert "primer3_explain" in primers
            if primers["pair_count"]:
                assert "pair_penalty" in primers["pairs"][0]
                assert "left_gc_percent" in primers["pairs"][0]
                assert "primer3_left_coordinate" in primers["pairs"][0]
                assert "primer3_right_coordinate" in primers["pairs"][0]

            pcrp = await call(
                s,
                "pcr_simulate",
                {"sequence_id": pcrt["sequence_id"], "forward_primer": fwd, "reverse_primer": rev, "anneal_min": 13},
            )
            assert pcrp["product_length_bp"] == 600
            assert "In-silico PCR" in pcrp["scope"]

            # Gibson assembly and primer design are separate: one simulates existing overlaps, one designs candidate tailed primers.
            gib_ids = []
            for i in (1, 2, 3):
                imp = await call(s, "sequence_import", {"path": str(ROOT / f"gibson{i}.fa"), "label": f"gib{i}"})
                gib_ids.append(imp["sequence_id"])
            gib = await call(s, "assembly_gibson", {"sequence_ids": gib_ids, "overlap_min": 25, "circular_only": True})
            assert gib["product_count"] >= 1

            design_ids = []
            for i in (1, 2, 3):
                imp = await call(s, "sequence_import", {"path": str(ROOT / f"design{i}.fa"), "label": f"design{i}"})
                design_ids.append(imp["sequence_id"])
            gd = await call(
                s,
                "gibson_primer_design",
                {"sequence_ids": design_ids, "overlap": 30, "target_tm": 60.0, "circular": True},
            )
            assert len(gd["designs"]) == 3
            assert "exact fragment/template" in gd["input_semantics"]
            assert len(gd["input_topologies"]) == 3
            assert all(d["forward_primer"]["full_length"] > d["forward_primer"]["annealing_length"] for d in gd["designs"])
            assert all(d["reverse_primer"]["full_length"] > d["reverse_primer"]["annealing_length"] for d in gd["designs"])
            assert all(d["simulated_pcr_product_length_bp"] > 500 for d in gd["designs"])

            # Golden Gate must reject ordinary restriction enzymes and accept BsaI geometry.
            assess = await call(s, "golden_gate_assess", {"sequence_ids": design_ids[:1], "enzymes": ["EcoRI", "BsaI"]})
            profiles = {x["name"]: x for x in assess["enzyme_profiles"]}
            assert profiles["EcoRI"]["golden_gate_type_iis_geometry"] is False
            assert profiles["BsaI"]["golden_gate_type_iis_geometry"] is True
            await call(
                s,
                "assembly_golden_gate",
                {"sequence_ids": design_ids[:2], "enzymes": ["EcoRI"], "circular_only": True},
                expect_error=True,
            )

            gg_ids = []
            for i in (1, 2, 3, 4):
                imp = await call(s, "sequence_import", {"path": str(ROOT / f"gg{i}.gb"), "label": f"gg{i}"})
                gg_ids.append(imp["sequence_id"])
            gg = await call(
                s,
                "assembly_golden_gate",
                {"sequence_ids": gg_ids, "enzymes": ["BsaI"], "circular_only": True, "max_products": 5},
            )
            assert gg["product_count"] == 1, gg

            status = await call(s, "tool_status", {})
            assert status["version"] == "1.2.0"

            print("PASS_EDIT_STALE_TRANSLATION_REMOVAL")
            print("PASS_CIRCULAR_RESTRICTION_ORIGIN")
            print("PASS_PRIMER3_QC")
            print("PASS_PCR", pcrp["product_length_bp"])
            print("PASS_GIBSON", gib["product_count"])
            print("PASS_GIBSON_PRIMER_DESIGN", len(gd["designs"]))
            print("PASS_GOLDEN_GATE_ASSESS")
            print("PASS_GOLDEN_GATE_BSAI", gg["product_count"])
            print("PASS_CONSTRUCT_VALIDATION")
            print("PASS_SPECIAL_CDS_QUALIFIER_HANDLING")
            print("PASS_FACTUAL_SCOPE_LABELS")
            print("ACCEPTANCE_V120_PASS")


asyncio.run(main())
