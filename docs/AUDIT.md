# 1.0 API / numerical audit

Completed as part of epic E10.

## Proposal §6 required features

All items in [`SUPPORT.md`](SUPPORT.md) “Required (1.0)” are implemented with tests
and, where applicable, SciPy fixtures listed in [`../fixtures/CATALOG.md`](../fixtures/CATALOG.md).

## Fixture catalog

The catalog maps every required numerical feature to a committed fixture or
law/behavioral suite. This is mapped fixture coverage, not a claim that every
public API path or scientific method is complete. The stale “gap” row for
lfilter/sosfilt/filtfilt was removed — those are covered by E4 fixtures.

End-to-end workflow: `smoke.e2e_butter_welch` (SOS → Welch).

## Benchmark receipts

Present under `benchmarks/receipts/` for Auto dispatch, FFT, FFT convolution,
direct operator, and SOS (`e10-sos.md`). Comparative SciPy timings:
[`../benchmarks/receipts/vs-scipy.md`](../benchmarks/receipts/vs-scipy.md)
(see [`PERFORMANCE.md`](PERFORMANCE.md)).

## Docs gate

- Support matrix: [SUPPORT.md](SUPPORT.md)
- Migration examples: [MIGRATION.md](MIGRATION.md)
- Intentional differences: [SEMANTICS.md](SEMANTICS.md)

## Ravel

`signal4s-ravel` provides `AlongAxis` / `ChannelBatch` with explicit
`ContiguousPolicy`. Core modules do not depend on ravel.
