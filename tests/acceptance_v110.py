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

ROOT = Path(r"C:\ProgramData\ChatGPT-Operator\workspace\gene_workbench_dist\acceptance_v110")
ROOT.mkdir(parents=True, exist_ok=True)
PYTHON = r"D:\bio-ai-workbench\envs\sequence\Scripts\python.exe"
SERVER = r"C:\ProgramData\ChatGPT-Operator\workspace\gene_workbench_dist\gene_workbench.py"


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


def prepare_inputs():
    cds = "ATG" + ("GCC" * 98) + "TAA"
    edit_seq = ("ACGT" * 25) + cds + ("TGCA" * 25)
    features = [
        SeqFeature(SimpleLocation(100, 400, strand=1), type="CDS", qualifiers={"label": ["good_cds"]}),
        SeqFeature(SimpleLocation(420, 480, strand=1), type="misc_feature", qualifiers={"label": ["tail_feature"]}),
    ]
    write_gb(ROOT / "edit_source.gb", edit_seq, False, features)

    pcr_seq = randseq(600, 7)
    write_gb(ROOT / "pcr_template.gb", pcr_seq, False, [])
    fwd = pcr_seq[:24]
    rev = str(Seq(pcr_seq[-24:]).reverse_complement())

    z = randseq(25, 101)
    x = randseq(25, 102)
    y = randseq(25, 103)
    a = randseq(70, 201)
    b = randseq(70, 202)
    c = randseq(70, 203)
    write_fa(ROOT / "gibson1.fa", z + a + x)
    write_fa(ROOT / "gibson2.fa", x + b + y)
    write_fa(ROOT / "gibson3.fa", y + c + z)

    write_gb(ROOT / "gg_backbone.gb", "cccGAATTCaaaGTCGACccc", True, [])
    write_gb(ROOT / "gg_insert.gb", "ggGAATTCaggtGTCGACgg", False, [])

    return fwd, rev


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


async def call(session, name, args):
    r = await session.call_tool(name, args)
    is_error = getattr(r, "isError", getattr(r, "is_error", False))
    if is_error:
        raise RuntimeError(f"{name} failed: {r}")
    return payload(r)


async def main():
    fwd, rev = prepare_inputs()
    params = StdioServerParameters(command=PYTHON, args=[SERVER])
    async with stdio_client(params) as (rs, ws):
        async with ClientSession(rs, ws) as s:
            await s.initialize()
            tools = await s.list_tools()
            names = [t.name for t in tools.tools]
            print("TOOL_COUNT", len(names))
            expected_new = {
                "sequence_insert", "sequence_delete", "sequence_diff", "feature_remap",
                "pcr_simulate", "assembly_gibson", "assembly_golden_gate",
                "assembly_ligation", "construct_validate"
            }
            missing = sorted(expected_new - set(names))
            assert not missing, missing

            src = await call(s, "sequence_import", {"path": str(ROOT / "edit_source.gb"), "label": "edit_src"})
            src_id = src["sequence_id"]
            val = await call(s, "construct_validate", {"sequence_id": src_id})
            assert val["overall_pass"] is True
            assert val["cds_count"] == 1
            assert val["cds_checks"][0]["frame_ok"] is True

            ins = await call(s, "sequence_insert", {
                "sequence_id": src_id, "position": 50, "insert_sequence": "GATTACAGAT", "label": "inserted"
            })
            assert ins["new_length_bp"] == ins["old_length_bp"] + 10
            ins_features = await call(s, "sequence_features", {"sequence_id": ins["sequence_id"]})
            cds_row = next(x for x in ins_features["features"] if x["type"] == "CDS")
            assert cds_row["start"] == 111 and cds_row["end"] == 410, cds_row

            dele = await call(s, "sequence_delete", {
                "sequence_id": ins["sequence_id"], "start": 50, "end": 59, "label": "restored"
            })
            assert dele["new_length_bp"] == dele["old_length_bp"] - 10

            diff = await call(s, "sequence_diff", {
                "source_sequence_id": src_id, "target_sequence_id": ins["sequence_id"]
            })
            assert diff["change_count"] >= 1
            assert diff["change_type_counts"]["insert"] >= 1

            remap = await call(s, "feature_remap", {
                "source_sequence_id": src_id, "target_sequence_id": ins["sequence_id"], "label": "remapped"
            })
            assert remap["features_mapped"] >= 2

            pcrt = await call(s, "sequence_import", {"path": str(ROOT / "pcr_template.gb"), "label": "pcr_template"})
            pcrp = await call(s, "pcr_simulate", {
                "sequence_id": pcrt["sequence_id"], "forward_primer": fwd, "reverse_primer": rev, "anneal_min": 13
            })
            assert pcrp["product_length_bp"] == 600, pcrp

            gib_ids = []
            for i in (1, 2, 3):
                imp = await call(s, "sequence_import", {"path": str(ROOT / f"gibson{i}.fa"), "label": f"gib{i}"})
                gib_ids.append(imp["sequence_id"])
            gib = await call(s, "assembly_gibson", {
                "sequence_ids": gib_ids, "overlap_min": 25, "circular_only": True, "max_products": 5
            })
            assert gib["product_count"] >= 1, gib
            assert any(p["topology"] == "circular" for p in gib["products"])

            bb = await call(s, "sequence_import", {"path": str(ROOT / "gg_backbone.gb"), "label": "ggbb"})
            ii = await call(s, "sequence_import", {"path": str(ROOT / "gg_insert.gb"), "label": "ggins"})
            gg_ids = [bb["sequence_id"], ii["sequence_id"]]
            gg = await call(s, "assembly_golden_gate", {
                "sequence_ids": gg_ids, "enzymes": ["EcoRI", "SalI"], "circular_only": True, "max_products": 5
            })
            assert gg["product_count"] >= 1, gg

            lig = await call(s, "assembly_ligation", {
                "sequence_ids": gg_ids, "enzymes": ["EcoRI", "SalI"], "circular_only": True, "max_products": 5
            })
            assert lig["product_count"] >= 1, lig

            product_id = gg["products"][0]["sequence_id"]
            product_val = await call(s, "construct_validate", {"sequence_id": product_id})
            assert product_val["overall_pass"] is True

            status = await call(s, "tool_status", {})
            assert status["version"] == "1.1.0"

            print("PASS_EDIT_REMAP")
            print("PASS_DIFF")
            print("PASS_PCR", pcrp["product_length_bp"])
            print("PASS_GIBSON", gib["product_count"])
            print("PASS_GOLDEN_GATE", gg["product_count"])
            print("PASS_LIGATION", lig["product_count"])
            print("PASS_VALIDATE")
            print("ACCEPTANCE_V110_PASS")


asyncio.run(main())
