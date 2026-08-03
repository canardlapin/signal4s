# Repository documentation

The public learning path is maintained in [`site-docs/`](../site-docs/README.md)
and is compiled by `sbt docs/tlSite`. This directory retains the detailed
reference and development record for the repository.

| Doc | Purpose |
|---|---|
| [SUPPORT.md](SUPPORT.md) | Detailed supported / deferred / platform matrix |
| [MIGRATION.md](MIGRATION.md) | Detailed SciPy / gsignal → signal4s maps and examples |
| [SEMANTICS.md](SEMANTICS.md) | Detailed intentional API/numeric differences |
| [PERFORMANCE.md](PERFORMANCE.md) | Scoped benchmark interpretation and receipts |
| [AUDIT.md](AUDIT.md) | Development audit; not part of the public site |
| [`../proposal.md`](../proposal.md) | Architecture and scope constitution |
| [`../fixtures/CATALOG.md`](../fixtures/CATALOG.md) | SciPy fixture coverage |

The site source deliberately excludes the audit, proposal, fixture generator,
and raw receipt material from its public navigation.
