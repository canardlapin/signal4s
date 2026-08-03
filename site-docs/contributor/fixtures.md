# Fixture-backed verification

Numeric parity is checked from committed fixtures generated with pinned
SciPy/NumPy versions. Tests do not invoke Python at runtime. The generator and
its environment are under `fixtures/generate/`; the catalog maps each fixture
to the implementation path and test suite.

For contributor work:

1. read [`fixtures/README.md`](https://github.com/canardlapin/signal4s/blob/main/fixtures/README.md),
2. inspect the mapping in [`fixtures/CATALOG.md`](https://github.com/canardlapin/signal4s/blob/main/fixtures/CATALOG.md),
3. regenerate only when the reference versions and tolerances are intentional,
4. run the JVM parity suite and the shared JVM/Scala.js suites.

Fixture parity is evidence for the declared cases. It is not a substitute for
the semantic laws, cross-platform tests, API documentation, or performance
measurements.
