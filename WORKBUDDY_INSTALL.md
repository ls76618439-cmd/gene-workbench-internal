# WorkBuddy installation protocol

This file is written for WorkBuddy. Follow it exactly. Do not redesign the installation.

## Goal

Install Gene Workbench locally on Windows, register the `gene-workbench` MCP in WorkBuddy, install the bundled Skill, and verify the local configuration.

## Preconditions

- Windows
- WorkBuddy already installed
- This GitHub repository is accessible to the user

Do not ask the user to install Python, edit JSON, configure ports, or manually copy Skill files.

## Install

From a clone of this repository, run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\bootstrap.ps1
```

The bootstrap will install GitHub CLI with winget if it is missing. If GitHub authentication is required for this private repository, ask the user to complete `gh auth login` once and then rerun the bootstrap.

Do not manually edit `~\.workbuddy\mcp.json` unless `bootstrap.ps1` reports a failure.

## Expected success

The script must print:

```text
GENE_WORKBENCH_INSTALL_OK
VERSION=1.2.0
```

Then verify:
- `%LOCALAPPDATA%\GeneWorkbench\GeneWorkbench-1.2.0.exe` exists
- `~\.workbuddy\mcp.json` contains `gene-workbench`
- `~\.workbuddy\skills\gene-workbench\SKILL.md` exists

If WorkBuddy was open during installation, restart or reload WorkBuddy once so it rediscovers the MCP tools.

After reload, call `tool_status` from `gene-workbench`. Expect version `1.2.0`. Then list tools and expect 22 tools. If both checks pass, report that installation and runtime discovery passed.

## Normal use

For DNA, plasmid, vector, primer, PCR, restriction digest, cloning, Gibson, Golden Gate, ligation, construct validation, or sequence-editing requests, use Gene Workbench when deterministic sequence computation would improve correctness.

Do not force Gene Workbench for purely conceptual questions. Do not force a cloning method when the user's goal and constraints do not determine one. Preserve tool warnings and distinguish computational results from experimental validation.
