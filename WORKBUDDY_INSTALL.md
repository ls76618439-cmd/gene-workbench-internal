# WorkBuddy installation protocol

This file is written for WorkBuddy. Follow it exactly. Do not redesign the installation.

## Goal

Install Gene Workbench locally on Windows, register the `gene-workbench` MCP in WorkBuddy, install the bundled Skill, and verify the local configuration.

## Preconditions

- Windows
- WorkBuddy already installed
- This GitHub repository is accessible to the user
- GitHub CLI (`gh`) is authenticated for this private repository

If GitHub authentication is the only missing precondition, ask the user to authenticate GitHub once, then continue. Do not ask the user to edit JSON or install Python.

## Install

From a clone of this repository, run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\bootstrap.ps1
```

Do not manually edit `~\.workbuddy\mcp.json` unless `bootstrap.ps1` reports a failure.

## Expected success

The script must print:

```text
GENE_WORKBENCH_INSTALL_OK
```

Then verify:
- `%LOCALAPPDATA%\GeneWorkbench\GeneWorkbench.exe` exists
- `~\.workbuddy\mcp.json` contains `gene-workbench`
- `~\.workbuddy\skills\gene-workbench\SKILL.md` exists

If WorkBuddy was open during installation, restart or reload WorkBuddy once so it rediscovers the MCP tools.

After reload, use `tool_status` from `gene-workbench`. If it succeeds, report that installation and runtime discovery both passed.

## Normal use

For DNA, plasmid, vector, primer, restriction digest, sequence edit or related sequence-analysis requests, use the `gene-workbench` MCP tools instead of manually reasoning over long raw sequences.
