---
name: gene-workbench
display_name: Gene Workbench
display_name_en: Gene Workbench
description: Analyze, edit, search, translate, digest, design primers for, and export DNA sequence files using the local gene-workbench MCP tools.
description_zh: 使用本地 Gene Workbench MCP 工具分析、编辑、检索、翻译、酶切、设计引物并导出 DNA 序列文件。
description_en: Analyze and edit DNA sequence files with local MCP tools.
version: 1.0.0
author: TREENEWBEE
---

# Gene Workbench

Use the `gene-workbench` MCP tools for DNA sequence work instead of manipulating long raw sequences in the model context.

## Preferred workflow

1. Import `.gb`, `.gbk`, `.genbank`, `.fa`, `.fasta`, `.fna`, `.fas`, or SnapGene `.dna` with `sequence_import`.
2. Keep the returned `sequence_id` and use it for later operations.
3. Use `sequence_info`, `sequence_features`, `sequence_find`, and `sequence_extract` for inspection.
4. Use `primer_design`, `restriction_analyze`, and `sequence_translate` for analysis.
5. Use `sequence_replace` for precise coordinate edits. Note that this V1 edit tool does not remap existing feature coordinates automatically.
6. Use `sequence_export` to produce GenBank or FASTA outputs.

For large sequences, do not paste the full sequence into chat unless necessary. Prefer IDs, coordinates, features, and short extracted regions.

Before a sequence-changing operation, summarize the intended coordinate range and replacement. After editing, verify length and relevant translated/restriction results before reporting success.
