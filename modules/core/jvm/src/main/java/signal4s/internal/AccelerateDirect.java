package signal4s.internal;

import java.lang.foreign.Arena;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SymbolLookup;
import java.lang.foreign.ValueLayout;
import java.lang.invoke.MethodHandle;
import java.nio.file.Path;

/**
 * macOS Accelerate implementation of zero-boundary full convolution.
 *
 * <p>vDSP_convD accepts one padded signal and performs the whole FIR in one
 * native call. The class is only selected after its macOS symbol lookup
 * succeeds; all other platforms retain the portable kernels.</p>
 */
final class AccelerateDirect {
  // vDSP_convD's documented vectorized path supports filter lengths <= 2044.
  static final int MAX_FILTER_LENGTH = 2044;

  private static volatile String diagnostic = "not initialized";
  private static volatile String dlopenDiagnostic = "not attempted";
  private static final MethodHandle VDSP_CONV_D = lookupConvD();

  /** Off-heap staging permits a regular FFM downcall without pinning heap arrays.
   * The buffers are reused per worker thread and replaced only when a new shape
   * exceeds their capacity.
   */
  private static final ThreadLocal<NativeScratch> SCRATCH =
      ThreadLocal.withInitial(NativeScratch::new);

  private AccelerateDirect() {}

  static boolean available() {
    return VDSP_CONV_D != null;
  }

  static String diagnostic() {
    return diagnostic;
  }

  static void convolveFull(double[] signal, double[] filter, double[] out) {
    final int n = signal.length;
    final int m = filter.length;
    if (m == 0 || out.length != n + m - 1 || m > MAX_FILTER_LENGTH) {
      throw new IllegalArgumentException("unsupported vDSP full-convolution shape");
    }

    final int leftPad = m - 1;
    final int paddedLength = n + 2 * leftPad;
    final NativeScratch scratch = SCRATCH.get();
    scratch.ensure(paddedLength, m, out.length);
    try {
      final long paddedBytes = (long) paddedLength * Double.BYTES;
      final long signalBytes = (long) n * Double.BYTES;
      final long filterBytes = (long) m * Double.BYTES;
      final long outBytes = (long) out.length * Double.BYTES;
      final MemorySegment padded = scratch.signal.asSlice(0, paddedBytes);
      final MemorySegment nativeFilter = scratch.filter.asSlice(0, filterBytes);
      final MemorySegment nativeOut = scratch.out.asSlice(0, outBytes);

      padded.fill((byte) 0);
      padded.asSlice((long) leftPad * Double.BYTES, signalBytes)
          .copyFrom(MemorySegment.ofArray(signal));
      nativeFilter.copyFrom(MemorySegment.ofArray(filter));

      // A positive filter stride computes correlation. Pointing to the final
      // tap with stride -1 computes y[t] = sum_j filter[j] * signal[t-j].
      VDSP_CONV_D.invokeExact(
          padded,
          1L,
          nativeFilter.asSlice((long) leftPad * Double.BYTES),
          -1L,
          nativeOut,
          1L,
          (long) out.length,
          (long) m);
      MemorySegment.ofArray(out).copyFrom(nativeOut);
    } catch (Throwable error) {
      throw new IllegalStateException("vDSP_convD failed", error);
    }
  }

  private static MethodHandle lookupConvD() {
    final String os = System.getProperty("os.name", "").toLowerCase();
    if (!os.contains("mac")) {
      diagnostic = "not macOS";
      return null;
    }
    Throwable loadFailure = null;
    try {
      // dyld can resolve framework names from its shared cache even though the
      // on-disk framework binary is a stub on recent macOS releases.
      System.loadLibrary("Accelerate");
      final MethodHandle resolved = resolve(SymbolLookup.loaderLookup());
      if (resolved != null) {
        diagnostic = "loaded via System.loadLibrary(Accelerate)";
        return resolved;
      }
      loadFailure = new UnsatisfiedLinkError("vDSP_convD absent from loader lookup");
    } catch (Throwable error) {
      loadFailure = error;
    }
    final MethodHandle viaDlopen = lookupViaDlopen(
        "/System/Library/Frameworks/Accelerate.framework/Accelerate");
    if (viaDlopen != null) {
      diagnostic = "loaded from dyld shared cache via dlopen";
      return viaDlopen;
    }
    try {
      final SymbolLookup accelerate =
          SymbolLookup.libraryLookup(
              Path.of("/System/Library/Frameworks/Accelerate.framework/Accelerate"),
              Arena.global());
      final MethodHandle resolved = resolve(accelerate);
      if (resolved != null) {
        diagnostic = "loaded via framework path";
        return resolved;
      }
      diagnostic =
          "vDSP_convD absent from framework lookup; dlopen: "
              + dlopenDiagnostic
              + "; prior load: "
              + loadFailure;
      return null;
    } catch (Throwable error) {
      diagnostic =
          "Accelerate load failed: "
              + error
              + "; dlopen: "
              + dlopenDiagnostic
              + "; prior load: "
              + loadFailure;
      return null;
    }
  }

  /** File-based libraryLookup rejects shared-cache-only framework binaries.
   * dlopen itself resolves their install names through dyld.
   */
  private static MethodHandle lookupViaDlopen(String path) {
    try {
      final Linker linker = Linker.nativeLinker();
      final SymbolLookup defaults = linker.defaultLookup();
      final MethodHandle dlopen = linker.downcallHandle(
          defaults.find("dlopen").orElseThrow(),
          FunctionDescriptor.of(
              ValueLayout.ADDRESS,
              ValueLayout.ADDRESS,
              ValueLayout.JAVA_INT));
      final MethodHandle dlsym = linker.downcallHandle(
          defaults.find("dlsym").orElseThrow(),
          FunctionDescriptor.of(
              ValueLayout.ADDRESS,
              ValueLayout.ADDRESS,
              ValueLayout.ADDRESS));

      final MemorySegment library = (MemorySegment) dlopen.invokeExact(utf8z(path), 1);
      if (library.address() == 0L) {
        return null;
      }
      final MemorySegment symbol =
          (MemorySegment) dlsym.invokeExact(library, utf8z("vDSP_convD"));
      if (symbol.address() == 0L) {
        return null;
      }
      return Linker.nativeLinker().downcallHandle(
          symbol,
          FunctionDescriptor.ofVoid(
              ValueLayout.ADDRESS,
              ValueLayout.JAVA_LONG,
              ValueLayout.ADDRESS,
              ValueLayout.JAVA_LONG,
              ValueLayout.ADDRESS,
              ValueLayout.JAVA_LONG,
              ValueLayout.JAVA_LONG,
              ValueLayout.JAVA_LONG));
    } catch (Throwable error) {
      dlopenDiagnostic = error.toString();
      return null;
    }
  }

  private static MemorySegment utf8z(String text) {
    // Downcalls require an off-heap pointer, not a heap segment. These two
    // process-lifetime strings are allocated once during static initialization.
    return Arena.global().allocateFrom(text);
  }

  private static MethodHandle resolve(SymbolLookup lookup) {
    final MemorySegment symbol = lookup.find("vDSP_convD").orElse(null);
    if (symbol == null) {
      return null;
    }
    return Linker.nativeLinker().downcallHandle(
        symbol,
        FunctionDescriptor.ofVoid(
            ValueLayout.ADDRESS,
            ValueLayout.JAVA_LONG,
            ValueLayout.ADDRESS,
            ValueLayout.JAVA_LONG,
            ValueLayout.ADDRESS,
            ValueLayout.JAVA_LONG,
            ValueLayout.JAVA_LONG,
            ValueLayout.JAVA_LONG));
  }

  private static final class NativeScratch {
    private Arena arena;
    private MemorySegment signal;
    private MemorySegment filter;
    private MemorySegment out;
    private int signalCapacity;
    private int filterCapacity;
    private int outCapacity;

    void ensure(int requiredSignal, int requiredFilter, int requiredOut) {
      if (arena != null
          && signalCapacity >= requiredSignal
          && filterCapacity >= requiredFilter
          && outCapacity >= requiredOut) {
        return;
      }
      if (arena != null) {
        arena.close();
      }
      arena = Arena.ofConfined();
      signalCapacity = requiredSignal;
      filterCapacity = requiredFilter;
      outCapacity = requiredOut;
      signal = arena.allocate((long) signalCapacity * Double.BYTES, Double.BYTES);
      filter = arena.allocate((long) filterCapacity * Double.BYTES, Double.BYTES);
      out = arena.allocate((long) outCapacity * Double.BYTES, Double.BYTES);
    }
  }
}
