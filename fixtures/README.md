# SciPy / gsignal fixtures

Versioned golden vectors for numeric evidence. **Scala tests load committed
fixture files only — they never invoke Python or R.**

## Layout

```text
fixtures/
  README.md
  schema/fixture.schema.json
  generate/                 # offline generator (not run in CI)
    requirements.txt        # pinned scipy/numpy
    generate_smoke.py
  data/                     # committed outputs
    smoke/
```

## Regenerating

Use a throwaway virtualenv. Do not commit the venv.

```bash
cd fixtures/generate
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python generate_smoke.py
```

Outputs land under `fixtures/data/`. Commit the JSON when intentional.

Pin bumps require regenerating fixtures and reviewing diffs. Record the new
`scipy` / `numpy` versions inside each fixture's `versions` object (also written
by the generator).

## Mapping notes (intentional API differences)

SciPy treats a signal as a bare array. `signal4s` carries axes on values:

| SciPy | signal4s mapping |
|---|---|
| `convolve(..., mode="full")` | `OutputRegion.Full` with causal or explicit `Kernel` origin |
| `convolve(..., mode="valid")` | `OutputRegion.Valid` |
| `convolve(..., mode="same")` | **Not** identical to a silent center. Map only via an explicit `Kernel` origin fixture that documents which SciPy alignment was used |
| bare coefficient vector | `Kernel(taps, zeroLagIndex)` |
| frequency as raw float | `Frequency` / `RadiansPerSample` / `CyclesPerSample` + `SampleRate` |

Fixtures that exercise `mode="same"` must include a `mapping` field describing
the origin convention used to compare against `signal4s`.

## CI policy

CI runs Scala tests against committed `fixtures/data/**`. It must not execute
`fixtures/generate/`.
