package signal4s.ravel

/** How axis vectors are presented to 1-D kernels.
  *
  * Contiguous copies are never implicit: callers choose a policy.
  */
enum ContiguousPolicy:
  /** Always gather into an owned contiguous [[gale.linalg.DVec]] before the op. */
  case CopyAlways
  /** Pass a strided view only when storage is already contiguous along the axis;
    * otherwise fail.
    */
  case RequireContiguous
