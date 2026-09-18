# Gene Workbench Internal

Internal Windows tool for using common molecular-biology sequence utilities from WorkBuddy through a local MCP server.

## What it does

Current V1 capabilities:
- Import GenBank / FASTA / SnapGene DNA files
- Sequence metadata and feature inspection
- Search and extract sequence regions
- Translation
- Restriction-site analysis
- Primer design
- Coordinate-based sequence replacement
- Export GenBank / FASTA

WorkBuddy uses the local `gene-workbench` MCP server over stdio. No server deployment is required.

## Easiest installation

Give WorkBuddy this repository URL and say:

> 请严格按照仓库根目录的 WORKBUDDY_INSTALL.md 安装 Gene Workbench，完成后验证 gene-workbench MCP 和 Skill 已安装成功。不要自行改安装方案。

For manual installation, download the latest `GeneWorkbench-Setup-*.exe` from Releases and run it.

## Repository layout

- `WORKBUDDY_INSTALL.md` - deterministic instructions for WorkBuddy
- `bootstrap.ps1` - one-command installer used by WorkBuddy
- `skill/gene-workbench/SKILL.md` - WorkBuddy routing skill
- `release-manifest.json` - expected release asset and checksum

## Internal distribution

This repository is private. A coworker needs read access to the repository and an authenticated GitHub CLI session for fully automatic installation.
