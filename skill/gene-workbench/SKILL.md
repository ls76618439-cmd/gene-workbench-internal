---
name: gene-workbench
display_name: Gene Workbench
display_name_en: Gene Workbench
description: Analyze, edit, search, translate, design primers, simulate PCR, analyze restriction sites, assemble DNA constructs, validate constructs, and export sequence files using the local gene-workbench MCP tools.
description_zh: 使用本地 Gene Workbench MCP 工具进行 DNA/质粒/载体结构分析、引物设计、PCR 模拟、酶切分析、序列编辑、Gibson/Golden Gate/连接组装、构建验证和序列导出。
description_en: Analyze and edit DNA sequence files and design common cloning workflows with local MCP tools.
version: 1.1.0
author: TREENEWBEE
---

# Gene Workbench

Use the `gene-workbench` MCP tools for DNA, plasmid, vector, primer, PCR, restriction digest, cloning, assembly, and sequence-editing tasks. Do not manually reason over long raw DNA sequences when a tool can do the operation deterministically.

## Routing

- Import: `sequence_import`
- Metadata / GC / topology: `sequence_info`
- Features / CDS / annotations: `sequence_features`
- Exact search / extraction: `sequence_find`, `sequence_extract`
- Translation: `sequence_translate`
- Primer design: `primer_design`
- PCR simulation: `pcr_simulate`
- Restriction sites: `restriction_analyze`
- Insert / delete / replace: `sequence_insert`, `sequence_delete`, `sequence_replace`
- Compare constructs: `sequence_diff`
- Annotation transfer after an externally-created sequence change: `feature_remap`
- Gibson assembly: `assembly_gibson`
- Golden Gate / type-IIS-style restriction assembly: `assembly_golden_gate`
- Restriction/ligation or direct ligation: `assembly_ligation`
- Final validation: `construct_validate`
- Export: `sequence_export`

## Preferred workflow

1. Import source sequence files and retain the returned `sequence_id`.
2. Inspect topology and relevant features before editing or assembly.
3. Use deterministic tools for primer, PCR, restriction, edit, and assembly calculations.
4. Treat every edit or assembly product as a new `sequence_id`; do not overwrite the source construct.
5. After an edit, use `sequence_diff` when the user needs a clear change summary.
6. After PCR or assembly, run `construct_validate` before declaring the construct ready.
7. Export the final construct as GenBank when annotations matter; use FASTA only when sequence alone is sufficient.

## Assembly notes

- `assembly_gibson` expects fragments that already contain compatible terminal overlaps.
- `assembly_golden_gate` and enzyme-assisted `assembly_ligation` take restriction enzyme names such as BsaI, BsmBI, EcoRI, or SalI.
- Multiple candidate products may be returned. Inspect topology, length, features, and validation results before selecting one.
- If a target was created outside the edit tools and its annotations need to follow sequence changes, use `feature_remap`.

For large sequences, keep raw DNA out of the chat when possible. Prefer sequence IDs, coordinates, features, short extracts, diffs, and exported files.
