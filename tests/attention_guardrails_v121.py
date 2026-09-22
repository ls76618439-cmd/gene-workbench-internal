from __future__ import annotations

import ast
from pathlib import Path

REPO = Path(r"D:\gene-workbench-internal")
SKILL = REPO / "skill" / "gene-workbench" / "SKILL.md"
SOURCE = REPO / "src" / "gene_workbench.py"

skill = SKILL.read_text(encoding="utf-8")
lines = skill.splitlines()
assert len(lines) <= 35, f"Skill too long: {len(lines)} lines"
assert len(skill) <= 3200, f"Skill too verbose: {len(skill)} chars"

# Detailed scientific limitations belong in tool descriptions/results, not repeated in the Skill.
for phrase in [
    "genome-wide",
    "hairpin",
    "type-IIS",
    "Primer3 core",
    "internal stop",
    "translation table",
]:
    assert phrase.lower() not in skill.lower(), f"Detailed limitation leaked into Skill: {phrase}"

assert "Do not force a tool call" in skill
assert "Do not force a cloning method" in skill
assert "tool's own `scope` and `warnings`" in skill

tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
tools: dict[str, str] = {}
for node in tree.body:
    if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        continue
    is_tool = False
    for dec in node.decorator_list:
        if isinstance(dec, ast.Call) and isinstance(dec.func, ast.Attribute):
            if isinstance(dec.func.value, ast.Name) and dec.func.value.id == "mcp" and dec.func.attr == "tool":
                is_tool = True
    if is_tool:
        tools[node.name] = ast.get_docstring(node) or ""

assert len(tools) == 22, f"Expected 22 MCP tools, got {len(tools)}"
for name, doc in tools.items():
    assert len(doc) >= 70, f"Tool description too thin: {name} ({len(doc)} chars)"

required_phrases = {
    "primer_design": ["Specificity outside the supplied template is not assessed"],
    "pcr_simulate": ["in-silico amplification", "not a claim"],
    "assembly_gibson": ["already contain compatible terminal overlaps", "does not design overlaps or primers"],
    "gibson_primer_design": ["supplied fragment order", "not assessed"],
    "golden_gate_assess": ["without selecting a preferred enzyme"],
    "assembly_golden_gate": ["compatible with this model", "caller chooses the enzyme"],
    "construct_validate": ["not experimental validation"],
}
for name, phrases in required_phrases.items():
    doc = tools[name].lower()
    for phrase in phrases:
        assert phrase.lower() in doc, f"{name} missing boundary phrase: {phrase}"

source = SOURCE.read_text(encoding="utf-8")
assert '"golden_gate_model_compatible"' in source
assert '"golden_gate_type_iis_geometry"' in source
assert '"computational_consistency_pass"' in source
assert '"overall_pass"' in source
assert 'deprecated compatibility alias' in source

print("SKILL_LINES", len(lines))
print("SKILL_CHARS", len(skill))
print("TOOL_COUNT", len(tools))
print("ATTENTION_GUARDRAILS_V121_PASS")
