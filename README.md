# Gene Workbench Internal

Internal Windows tool for common molecular-biology sequence and cloning operations from WorkBuddy through a local MCP server.

## Current V1.1 capabilities

- Import GenBank / FASTA / SnapGene DNA files
- Sequence metadata, GC, topology, feature and CDS inspection
- Exact sequence search and extraction
- Translation
- Primer3 primer design
- PCR simulation
- Restriction-site analysis
- Insert / delete / replace sequence edits with annotation remapping where possible
- Construct-to-construct diff
- Annotation remapping between related constructs
- Gibson assembly
- Golden Gate / restriction-and-ligation assembly
- Direct ligation assembly
- Construct validation including feature bounds and CDS reading-frame checks
- GenBank / FASTA export

WorkBuddy uses the local `gene-workbench` MCP server over stdio. No server deployment is required.

## Easiest installation

Give WorkBuddy this repository URL and say:

> 请严格按照仓库根目录的 WORKBUDDY_INSTALL.md 安装 Gene Workbench，完成后验证 gene-workbench MCP 和 Skill 已安装成功。不要自行改安装方案。

For manual installation, download the latest Company ZIP from Releases, extract it, and run `install.cmd`.

## Repository layout

- `WORKBUDDY_INSTALL.md` - deterministic instructions for WorkBuddy
- `bootstrap.ps1` - one-command installer used by WorkBuddy
- `skill/gene-workbench/SKILL.md` - WorkBuddy routing skill
- `src/gene_workbench.py` - MCP source
- `tests/acceptance_v110.py` - end-to-end acceptance
- `release-manifest.json` - pinned release asset and checksum

## Internal distribution

This repository is private. A coworker needs read access to the repository. The bootstrap installs GitHub CLI automatically when possible; the coworker only needs to complete GitHub authentication once if it is not already configured.
