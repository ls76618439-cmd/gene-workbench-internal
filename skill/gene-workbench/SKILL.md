---
name: gene-workbench
display_name: Gene Workbench
display_name_en: Gene Workbench
description: Use the local Gene Workbench MCP for sequence analysis, primer design, PCR simulation, restriction analysis, sequence editing, Gibson/Golden Gate/ligation assembly, construct consistency checks, and sequence export when deterministic sequence computation would improve correctness.
description_zh: 当确定性的序列计算能提高准确性时，使用本地 Gene Workbench MCP 进行序列分析、引物设计、PCR 模拟、酶切分析、序列编辑、Gibson/Golden Gate/连接组装、构建一致性检查和序列导出。
description_en: Use local deterministic sequence tools when they improve correctness for molecular-biology sequence work.
version: 1.2.0
author: TREENEWBEE
---

# Gene Workbench

Use Gene Workbench as a computation layer, not as a substitute for scientific judgment.

## Core principles

- Do not force a tool call for conceptual questions that can be answered accurately without sequence computation.
- Treat the tool map and usage patterns as navigation aids, not a mandatory workflow. A capable model may combine tools differently when the user's task and evidence support it.
- Do not force a cloning method. Respect an explicitly requested method. If the method is underdetermined, state the assumptions and present feasible options without calling one "best" unless the user supplied decision criteria.
- Treat tool output as computation on the supplied sequences and parameters. Do not describe computational feasibility as guaranteed wet-lab success.
- Preserve warnings, ambiguity, multiple candidate products, and scope limitations returned by tools.
- For long DNA sequences, prefer sequence IDs, coordinates, features, short extracts, diffs, and exported files over copying the full sequence into chat.

## Tool map

- Import: `sequence_import`
- Metadata / GC / topology: `sequence_info`
- Features / CDS / annotations: `sequence_features`
- Exact search / extraction: `sequence_find`, `sequence_extract`
- Translation: `sequence_translate`
- Primer candidates: `primer_design`
- PCR simulation: `pcr_simulate`
- Restriction sites: `restriction_analyze`
- Insert / delete / replace: `sequence_insert`, `sequence_delete`, `sequence_replace`
- Compare constructs: `sequence_diff`
- Annotation transfer: `feature_remap`
- Gibson assembly from fragments that already contain overlaps: `assembly_gibson`
- Candidate overlap-tailed Gibson PCR primers for a specified fragment order: `gibson_primer_design`
- Golden Gate enzyme/site facts without choosing an enzyme: `golden_gate_assess`
- Golden Gate assembly with compatible type-IIS cut geometry: `assembly_golden_gate`
- Ordinary restriction/ligation or direct ligation: `assembly_ligation`
- Computational construct consistency checks: `construct_validate`
- Export: `sequence_export`

## Scientific interpretation

- `primer_design` returns Primer3 candidates and quality metrics for the supplied template. Do not claim genome-wide or host-wide specificity unless an external reference/off-target search was actually performed.
- `pcr_simulate` models annealing and product formation on the supplied template; it does not establish real PCR efficiency.
- `gibson_primer_design` designs overlap-tailed primers for the supplied fragment order and simulates those primer pairs on the supplied templates. It does not prove experimental assembly efficiency.
- `golden_gate_assess` reports recognition/cut geometry and site counts. It intentionally does not select a preferred enzyme.
- `assembly_golden_gate` is for enzymes whose cut geometry is compatible with Golden Gate-style type-IIS assembly. Use `assembly_ligation` for ordinary restriction/ligation.
- `construct_validate` checks sequence/annotation/CDS consistency. A passing result is not evidence that a physical construct is correct, expresses properly, or functions biologically.
- After sequence edits, note any dropped/remapped features or removed stale translation qualifiers. Recompute or review annotations that intersect the edited region.

## Reliable usage patterns

- Inspect topology and relevant annotations before coordinate-sensitive edits or restriction analysis.
- Keep each edit or assembly product as a new `sequence_id`; preserve the source construct.
- Use `sequence_diff` when the user needs an explicit record of changes.
- When an assembly returns multiple products, report the candidates and discriminating facts rather than silently choosing one.
- Before exporting a final computational construct, use `construct_validate` when CDS or annotation integrity matters.
