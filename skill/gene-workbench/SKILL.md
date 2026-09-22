---
name: gene-workbench
display_name: Gene Workbench
display_name_en: Gene Workbench
description: Use the local Gene Workbench MCP when deterministic DNA sequence computation can improve correctness for sequence analysis, editing, primer/PCR work, cloning assembly, construct consistency checks, or export.
description_zh: 当确定性的 DNA 序列计算能提高准确性时，使用本地 Gene Workbench MCP 进行序列分析、编辑、引物/PCR、克隆组装、构建一致性检查或导出。
description_en: Use local deterministic DNA sequence tools when they improve correctness.
version: 1.2.1
author: TREENEWBEE
---

# Gene Workbench

Use Gene Workbench as a computation layer, not as a substitute for scientific judgment.

- Do not force a tool call for conceptual questions that do not need sequence computation.
- Do not force a cloning method. Respect the user's requested method; when the method is underdetermined, keep alternatives open and state assumptions.
- Treat tool output as computation on supplied sequences and parameters, not as proof of wet-lab success.
- Preserve warnings, ambiguity, multiple candidate products, and each tool's returned `scope`.
- Prefer sequence IDs, coordinates, features, short extracts, diffs, and exported files over pasting long DNA sequences.

## Capability groups

- **Analysis:** sequence metadata, features, search, extraction, translation.
- **Editing:** insert, delete, replace, diff, annotation remap.
- **Primer / PCR:** primer candidate design and in-silico PCR.
- **Cloning:** restriction analysis, Gibson design/assembly, Golden Gate assessment/assembly, ligation.
- **Validation / output:** computational construct consistency checks and sequence export.

Choose tools from their names and descriptions. Treat the tool's own `scope` and `warnings` as the source of detailed limitations.
