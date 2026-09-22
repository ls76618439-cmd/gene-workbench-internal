# Gene Workbench Internal

Internal Windows tool for molecular-biology sequence analysis and common cloning computations from WorkBuddy through a local MCP server.

## Current V1.2.1 capabilities

- Import GenBank / FASTA / SnapGene DNA files
- Sequence metadata, GC, topology, feature and CDS inspection
- Exact sequence search, extraction, and translation
- Primer3 primer design with quality metrics via the bundled `primer3_core.exe` executable
- PCR simulation on the supplied template
- Restriction-site analysis with circular-topology handling
- Insert / delete / replace sequence edits with annotation remapping where possible
- Construct-to-construct diff
- Annotation remapping between related constructs, with stale sequence-derived qualifiers marked/removed
- Gibson assembly for fragments that already contain terminal overlaps
- Candidate overlap-tailed PCR primer design for Gibson assembly in a specified fragment order
- Golden Gate enzyme/site assessment without choosing a preferred enzyme
- Golden Gate assembly restricted to compatible type-IIS cut geometry
- Ordinary restriction/ligation and direct ligation assembly
- Computational construct validation including feature bounds, CDS frame, internal stops, codon-table checks, and translation-annotation consistency
- GenBank / FASTA export

WorkBuddy uses the local `gene-workbench` MCP server over stdio. No server deployment is required.

## Scientific-use principle

Gene Workbench separates deterministic computation from scientific judgment. Tools return calculations, candidates, warnings, and scope limits; the Skill does not force a cloning method or present computational feasibility as guaranteed wet-lab success.

## Easiest installation

Give WorkBuddy this repository URL and say:

> 閻犲洩娓瑰妤呭冀閸忕厧鐦婚柣鎾楀倻娉㈤幖瀛樻尰閻楁挳鎯勯鑲╃Э闁?WORKBUDDY_INSTALL.md 閻庣懓顦抽ˉ?Gene Workbench闁挎稑鑻悾顒勫箣閹邦剚鍊靛Δ鐘茬焷閻?gene-workbench MCP 闁?Skill 鐎瑰憡褰冮悾銊ф啑閸涱喖鐏囬柛鏃傚枂閳ь剙鍊风粭澶屾啺娴ｈ棄娈伴悶娑樻湰閺佽偐鈧懓顦抽ˉ濠囧棘鐟欏嫷鏀抽柕?

For manual installation, download the latest Company ZIP from Releases, extract it, and run `install.cmd`.

## Repository layout

- `WORKBUDDY_INSTALL.md` - deterministic instructions for WorkBuddy
- `bootstrap.ps1` - one-command installer used by WorkBuddy
- `skill/gene-workbench/SKILL.md` - WorkBuddy routing and scientific-interpretation guidance
- `src/gene_workbench.py` - MCP source
- `tests/acceptance_v121.py` - v1.2 end-to-end scientific acceptance
- `release-manifest.json` - pinned release asset and checksum
- `packaging/` - installer helper scripts

## Internal distribution

This repository is private. A coworker needs read access to the repository. The bootstrap installs GitHub CLI automatically when possible; the coworker only needs to complete GitHub authentication once if it is not already configured.
