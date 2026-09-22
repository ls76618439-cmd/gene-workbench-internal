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

> 闁荤姴娲╁〒鐟邦嚕濡ゅ懎鍐€闁稿繒鍘ч惁濠氭煟閹炬鍊诲▔銏ゅ箹鐎涙ɑ灏伴柣妤佹尦閹嫰顢欓懖鈺冃梺?WORKBUDDY_INSTALL.md 闁诲海鎳撻ˇ鎶剿?Gene Workbench闂佹寧绋戦懟顖炴偩椤掑嫬绠ｉ柟閭﹀墯閸婇潧螖閻樿尙鐒烽柣?gene-workbench MCP 闂?Skill 閻庣懓鎲¤ぐ鍐偩閵娧勫晳闁告侗鍠栭悘鍥煕閺冨倸鏋傞柍褜鍓欓崐椋庣箔婢跺本鍟哄ù锝堟濞堜即鎮跺☉妯绘拱闁轰浇鍋愰埀顒傛嚀椤︽娊藟婵犲洤妫橀悷娆忓閺€鎶芥煏?

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

This repository is public so coworkers can install without GitHub authentication. The bootstrap downloads the pinned public Release asset directly, verifies SHA256, installs the local MCP runtime, and registers the bundled Skill with WorkBuddy.
