package signal4s.internal;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;

/**
 * JVM Vector-API dense FIR scatter: {@code out[i+j] += x[i] * h[j]}.
 *
 * Enabled only when the preferred species is wide enough that we are not on a
 * software-emulated path (Apple Silicon reports SPECIES_256 but emulates it
 * slowly; preferred width there is 2).
 */
final class VectorScatter {
  private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;
  private static final boolean ENABLE = computeEnable();

  private VectorScatter() {}

  private static boolean computeEnable() {
    String arch = System.getProperty("os.arch", "").toLowerCase();
    // aarch64 preferred length is 2; SPECIES_256 is emulated and much slower.
    if (arch.contains("aarch64") || arch.contains("arm64") || arch.equals("arm")) {
      return false;
    }
    return SPECIES.length() >= 4;
  }

  static void scatterFull(double[] x, double[] h, double[] out) {
    final int n = x.length;
    final int m = h.length;
    final int speciesLen = SPECIES.length();
    final int mBound = SPECIES.loopBound(m);
    for (int i = 0; i < n; i++) {
      final double xi = x[i];
      if (xi == 0.0) continue;
      final DoubleVector xv = DoubleVector.broadcast(SPECIES, xi);
      int j = 0;
      for (; j < mBound; j += speciesLen) {
        DoubleVector hv = DoubleVector.fromArray(SPECIES, h, j);
        DoubleVector ov = DoubleVector.fromArray(SPECIES, out, i + j);
        ov.add(hv.mul(xv)).intoArray(out, i + j);
      }
      for (; j < m; j++) {
        out[i + j] += xi * h[j];
      }
    }
  }

  static int speciesLength() {
    return SPECIES.length();
  }

  static boolean worthwhile() {
    return ENABLE;
  }
}
