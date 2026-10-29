#include "copyright.i"

#pragma once

//=============================================================================================
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// KAMD (Kinetically Accelerated Molecular Dynamics) GPU Kernels
//
// Authors: Zhen Huang, Haixin Wei, Ray Luo
//          Department of Molecular Biology and Biochemistry, UC Irvine
//
//=============================================================================================
//
// Overview
// --------
// KAMD accelerates biomolecular sampling by periodically resetting first-shell solvent
// velocities in a physically motivated manner (pseudo-Brownian dynamics).  The original
// CPU implementation downloaded all coordinates, velocities, and forces from the GPU,
// performed first-shell detection via brute-force O(N_ref * N_solvent) loops, reassigned
// velocities on the host, and uploaded everything back — transferring ~2.8 MB per reset
// event and incurring 66% overhead.
//
// This GPU implementation eliminates nearly all host-device data transfer by:
//   (a) reusing PMEMD's existing cell decomposition (neighbor list infrastructure) for
//       spatial proximity search entirely on-device, and
//   (b) modifying velocity buffers in-place on the GPU.
//
// The only KAMD-specific host/device data transfer is:
//   - Momentum monitoring: currently a fixed 18 doubles downloaded per check
//     (144 bytes; first six groups in the current monitor buffer layout)
//   - Velocity redirection: ngroups*3 doubles uploaded per reset event
//     (24*ngroups bytes; e.g., 96 bytes for 4 groups, 144 bytes for 6 groups)
// This is still nearly 4 orders of magnitude smaller than the original
// CPU-side pseudoBD path (~2.8 MB per reset event).
// Measured overhead: < 0.5% vs. standard GPU MD (compared to 66% for the CPU path).
//
//
// How PMEMD's Neighbor List (NL) Cell Decomposition Is Used
// ---------------------------------------------------------
// PMEMD partitions the simulation box into a 3D grid of cells for its PME nonbond
// calculation.  Each cell is approximately (cutoff + skin) Angstroms wide — typically
// ~10 A for an 8 A nonbond cutoff with 2 A skin.  KAMD piggybacks on this existing
// spatial decomposition rather than building its own data structure.
//
// Key NL data structures accessed by the KAMD kernels (all in device memory via cSim):
//
//   pImageAtom[sorted_idx]        -> original (topology-order) atom index.
//                                    The NL engine sorts atoms into cells; this array
//                                    maps a sorted position back to the original index.
//
//   pImageAtomLookup[orig_idx]    -> sorted position in the image arrays.
//                                    Inverse of pImageAtom.  Given an original atom index
//                                    (e.g. from the reference group definition), this
//                                    gives its current position in the spatially-sorted
//                                    arrays so we can read its coordinates, velocity, etc.
//
//   pImageCellID[sorted_idx]      -> packed cell coordinates (cx, cy, cz) for the atom.
//                                    Encoded as a single unsigned int using bit fields:
//                                      cx = cellID & CELL_ID_MASK          (bits 0-5)
//                                      cy = (cellID >> CELL_IDY_SHIFT) & CELL_ID_MASK
//                                      cz = (cellID >> CELL_IDZ_SHIFT) & CELL_ID_MASK
//                                    CELL_ID_MASK = 0x3f (6 bits per dimension, max 64 cells).
//
//   pNLNonbondCellStartEnd[linearCell]  -> uint2{start, end} atom range in sorted order.
//                                    All atoms in cell linearCell occupy the contiguous
//                                    range [start, end) in the image arrays.  The linear
//                                    cell index is computed as:
//                                      linearCell = cz * ycells * xcells + cy * xcells + cx
//
//   xcells, ycells, zcells        -> number of cells along each box dimension.
//
//   pImageX/Y/Z[sorted_idx]       -> full-precision (double) Cartesian coordinates in the
//                                    image (sorted) order.  These are absolute coordinates
//                                    suitable for direct distance computation with minimum-
//                                    image convention.  Note: pAtomXYSP is cell-relative
//                                    and NOT usable for cross-cell distance calculations.
//
//   pImageVelX/Y/Z[sorted_idx]    -> velocities in image order (half-step for leapfrog).
//   pImageMass[sorted_idx]        -> atomic masses in image order.
//   pImageInvMass[sorted_idx]     -> inverse masses (1/m) in image order.
//
//   pForceXAccumulator[sorted_idx]   -> fixed-point (long long int) bonded force accumulator.
//   pNBForceXAccumulator[sorted_idx] -> fixed-point nonbond force accumulator (separate in NTP).
//                                    Forces are stored as PMEAccumulator (long long int) and
//                                    must be descaled by ONEOVERFORCESCALE to get real units.
//                                    In NTP simulations (non-MPI), bonded and nonbond forces
//                                    are stored in separate accumulators; in NVT, both point
//                                    to the same memory, so the splitForces flag controls
//                                    whether the nonbond accumulator is added.
//
// Cell neighbor search strategy:
//   Since KAMD's first-shell cutoff (3 A for surface marking, 5 A for group assignment) is
//   smaller than the cell width (~10 A), searching the home cell plus its 26 immediate
//   neighbors (a 3x3x3 cube) guarantees that ALL atoms within 5 A are found.  This avoids
//   the need for a dedicated spatial index — we simply iterate over the 27 cells surrounding
//   each reference atom, and for each cell, loop over the contiguous atom range given by
//   pNLNonbondCellStartEnd.  Periodic boundary wrapping of cell indices handles atoms near
//   box edges.
//
//
// Kernels
// -------
// 1. kKAMD_FindFirstShell_kernel  — Spatial search using NL cells.  One thread per reference
//    atom.  Searches 27-cell neighborhood.  Marks surface reference atoms (pRecord == -1)
//    and stores closest-group assignments race-safely in pPacked via atomicMin on a 64-bit
//    word that bundles (dist^2 bits, group).
//
// 2. kKAMD_ExtractGroup_kernel    — Post-pass: copy the group from pPacked into pRecord
//    (preserving the -1 surface marks).
//
// 3. kKAMD_ComputeMomenta_kernel  — GPU reduction of per-group translational momentum.
//    One thread per reference atom.  The current monitor path clears/downloads the
//    first 18 doubles (6 groups * 3 xyz) from the device momenta buffer.
//
// 4. kKAMD_RedirectVelocities_kernel — In-place velocity redirection on GPU.  One thread
//    per atom.  Preserves atomic speed while changing direction.  Includes leapfrog half-step
//    force correction: advance v to full-step, redirect, step back to half-step.
//
// 5. kKAMD_ClearBuffers_kernel    — Reset pRecord to 0 and pPacked to (bits(cut^2) << 32).
// 6. kKAMD_ClearMomenta_kernel    — Zero the 18-double momenta buffer.
//
// These kernels replace the expensive CPU-side first-shell search and velocity
// redirection path that previously required full coordinate/velocity/force transfers.
//=============================================================================================

//---------------------------------------------------------------------------------------------
// 27-cell offset table for the 3D neighborhood search.
//
// For a reference atom in cell (cx, cy, cz), we must search all cells within one cell-width
// in each dimension to guarantee coverage of the 5 A KAMD cutoff (cell width ~10 A > 5 A).
// This is the standard 3x3x3 stencil: the home cell (0,0,0) plus 26 neighbors.
// Periodic boundary conditions are applied when cell indices go below 0 or above xcells-1.
//---------------------------------------------------------------------------------------------
__constant__ int kKAMD_CellOffset[27][3] = {
  { 0,  0,  0},
  { 1,  0,  0}, {-1,  0,  0},
  { 0,  1,  0}, { 0, -1,  0},
  { 0,  0,  1}, { 0,  0, -1},
  { 1,  1,  0}, {-1,  1,  0}, { 1, -1,  0}, {-1, -1,  0},
  { 1,  0,  1}, {-1,  0,  1}, { 1,  0, -1}, {-1,  0, -1},
  { 0,  1,  1}, { 0, -1,  1}, { 0,  1, -1}, { 0, -1, -1},
  { 1,  1,  1}, {-1,  1,  1}, { 1, -1,  1}, {-1, -1,  1},
  { 1,  1, -1}, {-1,  1, -1}, { 1, -1, -1}, {-1, -1, -1}
};

//---------------------------------------------------------------------------------------------
// kKAMD_FindFirstShell_kernel
//
// Purpose:
//   Identify first-shell solvent molecules around each reference group and mark surface
//   (solvent-exposed) reference-group atoms.  This is the spatial proximity search that
//   determines WHICH atoms will have their velocities redirected.
//
// How it uses the NL cell decomposition:
//   Each thread is assigned one reference-group atom.  Instead of computing distances to all
//   N_solvent atoms (O(N_ref * N_solvent) brute force), we exploit PMEMD's cell decomposition:
//
//   1. Look up the reference atom's sorted position via pImageAtomLookup[orig_r].
//   2. Read that atom's cell coordinates from pImageCellID[sorted_r], unpacking the packed
//      cell indices (cx, cy, cz) via bit masking (CELL_ID_MASK, CELL_IDY/Z_SHIFT).
//   3. Iterate over all 27 cells in the 3x3x3 neighborhood (home + 26 neighbors), applying
//      periodic wrapping when cell indices fall outside [0, xcells/ycells/zcells).
//   4. For each neighbor cell, compute the linear cell index and use
//      pNLNonbondCellStartEnd[linearCell] to get the contiguous range [start, end) of atoms
//      in that cell (atoms are sorted by cell during neighbor list construction).
//   5. Loop over atoms in that range, using pImageAtom[s] to recover the original atom index
//      and filtering to the active bath range for the region type
//      (solution bath or transmembrane bath).
//   6. Compute the distance using full-precision Cartesian coordinates pImageX/Y/Z (NOT the
//      cell-relative pAtomXYSP, which cannot be used for cross-cell distances) with the
//      minimum-image convention.  For orthogonal boxes this is a fast box-subtract; for
//      non-orthogonal boxes (e.g. truncated octahedron) we go through fractional coordinates
//      via the recip/ucell matrices, matching the convention in kReImageCoord.h.
//
//   This reduces the work per reference atom from O(N_solvent) to O(N_atoms_in_27_cells),
//   which is typically a few hundred atoms — a massive reduction for large systems.
//
// Two distance thresholds:
//   - r < 3 A (r2 < 9.0):   mark the reference atom as "surface" (pRecord[orig_r] = -1).
//                             Surface ref atoms will have their velocities redirected too.
//   - r < 5 A (r2 < kamd_cut2 = 25.0):  assign the solvent atom to its closest reference
//                             group (recorded in pPacked, extracted to pRecord post-pass).
//
// Race condition handling:
//   Multiple reference atoms may find the same solvent atom within cutoff and race to claim
//   it.  We pack {dist^2, groupNum} into a single 64-bit word — high 32 bits hold the IEEE
//   bit pattern of the float distance, low 32 bits hold the group number.  Because positive
//   floats compare correctly as unsigned ints, atomicMin on the packed word atomically keeps
//   the smallest distance AND its group together, eliminating the (CAS-dist, exch-group) race
//   pattern entirely.  After FindFirstShell completes, kKAMD_ExtractGroup_kernel copies the
//   group from pPacked into pRecord (preserving the -1 surface marks written separately).
//
// Thread mapping:
//   One thread per reference-group atom (totalRefAtoms threads total).
//   Block size: 256 threads, 1 block per SM max occupancy.
//
// Arguments:
//   pRecord         - per-atom flag (0=unaffected, -1=surface ref).  Solvent group numbers
//                     get extracted into pRecord by kKAMD_ExtractGroup_kernel after this kernel.
//   pPacked         - work array: (float-bits(dist^2) << 32) | groupNum, atomicMin'd per atom.
//                     Init to (float-bits(kamd_cut2) << 32) | 0.
//   totalRefAtoms   - number of reference-group atoms (sum over all active groups)
//   pRefAtomInfo    - per-ref-atom: int2{original_0based_index, group_number_1toN}
//   isNTP           - 1 if constant-pressure (read box from NTPData), 0 if constant-volume
//
// Box-shape handling:
//   The min-image convention uses the recip and ucell matrices.  For orthogonal boxes only the
//   diagonal entries are nonzero — that path collapses to a per-axis box subtraction.  For
//   non-orthogonal boxes (e.g. truncated octahedron / triclinic) the off-diagonal entries are
//   needed and we go through fractional coordinates: df = recip * d, df -= round(df),
//   d_minimg = ucell * df.  cSim.is_orthog selects the path; same convention as kReImageCoord.h
//   and kNL_PencilHash.h.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(256, 1)
kKAMD_FindFirstShell_kernel(
    int*           pRecord,
    unsigned long long* pPacked,
    int            totalRefAtoms,
    int2*          pRefAtomInfo,
    int            isNTP)
{
  // Constant-per-run KAMD parameters from cSim:
  //   bath ranges (0-based, inclusive) for solution and transmembrane regions,
  //   plus the squared distance cutoffs that select shell membership and
  //   surface marking.  Membrane regions use a single membrane_cut for both.
  const int   sol_start  = cSim.KAMDSolventStart;
  const int   sol_end    = cSim.KAMDSolventEnd;
  const int   mem_start  = cSim.KAMDMembraneStart;
  const int   mem_end    = cSim.KAMDMembraneEnd;
  const float shell_cut2   = cSim.KAMDShellCut2;
  const float surface_cut2 = cSim.KAMDSurfaceCut2;
  const float mem_cut2     = cSim.KAMDMembraneCut2;
  // Read box matrices for minimum-image convention.
  //   ucell is upper-triangular (box vectors in row-major: ucell[row][col]).
  //   recip is lower-triangular (inverse-transpose of ucell, same indexing).
  // For NTP the box changes each step → read from pNTPData (flat row-major at i*3+j).
  // For NVT/NVE the box is constant → use cSim.ucell/recip.
  // We hoist the 12 nonzero entries (6 of ucell + 6 of recip) once per thread.
  double u00, u01, u02, u11, u12, u22;
  double r00, r10, r11, r20, r21, r22;
  if (isNTP) {
    u00 = cSim.pNTPData->ucell[0]; u01 = cSim.pNTPData->ucell[1]; u02 = cSim.pNTPData->ucell[2];
    u11 = cSim.pNTPData->ucell[4]; u12 = cSim.pNTPData->ucell[5];
    u22 = cSim.pNTPData->ucell[8];
    r00 = cSim.pNTPData->recip[0];
    r10 = cSim.pNTPData->recip[3]; r11 = cSim.pNTPData->recip[4];
    r20 = cSim.pNTPData->recip[6]; r21 = cSim.pNTPData->recip[7]; r22 = cSim.pNTPData->recip[8];
  } else {
    u00 = cSim.ucell[0][0]; u01 = cSim.ucell[0][1]; u02 = cSim.ucell[0][2];
    u11 = cSim.ucell[1][1]; u12 = cSim.ucell[1][2];
    u22 = cSim.ucell[2][2];
    r00 = cSim.recip[0][0];
    r10 = cSim.recip[1][0]; r11 = cSim.recip[1][1];
    r20 = cSim.recip[2][0]; r21 = cSim.recip[2][1]; r22 = cSim.recip[2][2];
  }
  bool is_orthog = cSim.is_orthog;

  // Each thread processes one reference-group atom
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= (unsigned int)totalRefAtoms) return;

  // Read this reference atom's original index and group assignment
  int2 refInfo    = pRefAtomInfo[tid];
  int orig_r      = refInfo.x;        // original (topology-order) atom index, 0-based
  int groupNum    = refInfo.y;        // reference group number, 1..ngroups

  // Per-region-type bath range and cutoffs.
  //   region_type == KAMD_RT_SOLUTION:     bath = solvent range,  shell_cut for assignment, surface_cut for surface
  //   region_type == KAMD_RT_TRANSMEMBRANE: bath = membrane range, single membrane_cut for both
  int   regionType = cSim.pKAMDRegionType[groupNum - 1];
  int   bath_start;
  int   bath_end;
  float assign_cut2;
  float surf_cut2;
  if (regionType == 2) {  // KAMD_RT_TRANSMEMBRANE
    bath_start  = mem_start;
    bath_end    = mem_end;
    assign_cut2 = mem_cut2;
    surf_cut2   = mem_cut2;
  } else {                // KAMD_RT_SOLUTION
    bath_start  = sol_start;
    bath_end    = sol_end;
    assign_cut2 = shell_cut2;
    surf_cut2   = surface_cut2;
  }

  // --- NL lookup step 1: original index -> sorted position ---
  // pImageAtomLookup is the inverse of pImageAtom: given an original atom index, it returns
  // the atom's current position in the spatially-sorted image arrays.
  unsigned int sorted_r = cSim.pImageAtomLookup[orig_r];

  // Read full-precision Cartesian coordinates of the reference atom.
  // IMPORTANT: We use pImageX/Y/Z (absolute Cartesian, double precision), NOT pAtomXYSP
  // (single-precision, cell-relative offsets).  pAtomXYSP stores coordinates relative to
  // the cell origin, which cannot be used for distance calculations across cell boundaries.
  double xr = cSim.pImageX[sorted_r];
  double yr = cSim.pImageY[sorted_r];
  double zr = cSim.pImageZ[sorted_r];

  // --- NL lookup step 2: extract cell coordinates from packed cell ID ---
  // pImageCellID stores (cx, cy, cz) packed into a single unsigned int using 6 bits each:
  //   cx = bits [0:5],  cy = bits [6:11],  cz = bits [12:17]
  unsigned int cellID_r = cSim.pImageCellID[sorted_r];
  int cx_r = cellID_r & CELL_ID_MASK;
  int cy_r = (cellID_r >> CELL_IDY_SHIFT) & CELL_ID_MASK;
  int cz_r = (cellID_r >> CELL_IDZ_SHIFT) & CELL_ID_MASK;

  // --- NL lookup step 3: iterate over 27-cell neighborhood ---
  // The 3x3x3 stencil covers all cells within 1 cell-width (~10 A) of the home cell.
  // Since our cutoff (5 A) < cell width, this guarantees all atoms within cutoff are found.
  for (int nc = 0; nc < 27; nc++) {
    int cx_n = cx_r + kKAMD_CellOffset[nc][0];
    int cy_n = cy_r + kKAMD_CellOffset[nc][1];
    int cz_n = cz_r + kKAMD_CellOffset[nc][2];

    // Periodic boundary wrapping: if cell index goes below 0 or above max, wrap around
    if (cx_n < 0) cx_n += cSim.xcells;
    if (cx_n >= cSim.xcells) cx_n -= cSim.xcells;
    if (cy_n < 0) cy_n += cSim.ycells;
    if (cy_n >= cSim.ycells) cy_n -= cSim.ycells;
    if (cz_n < 0) cz_n += cSim.zcells;
    if (cz_n >= cSim.zcells) cz_n -= cSim.zcells;

    // --- NL lookup step 4: linear cell index -> atom range ---
    // pNLNonbondCellStartEnd[linearCell] gives uint2{start, end}: all atoms in this cell
    // occupy sorted positions [start, end) in the image arrays.
    int linearCell = cz_n * cSim.ycells * cSim.xcells +
                     cy_n * cSim.xcells + cx_n;
    uint2 cellRange = cSim.pNLNonbondCellStartEnd[linearCell];

    // --- NL lookup step 5: iterate over atoms in this cell ---
    for (unsigned int s = cellRange.x; s < cellRange.y; s++) {
      // pImageAtom[s] maps sorted position back to original atom index
      unsigned int orig_s = cSim.pImageAtom[s];

      // Filter: only process bath atoms for this group's region type
      // (solvent for solution regions, membrane for transmembrane regions).
      if ((int)orig_s < bath_start || (int)orig_s > bath_end) {
        continue;
      }

      // Read solvent atom coordinates (full-precision Cartesian, image order)
      double xs = cSim.pImageX[s];
      double ys = cSim.pImageY[s];
      double zs = cSim.pImageZ[s];

      double dx = xr - xs;
      double dy = yr - ys;
      double dz = zr - zs;

      // Minimum-image convention.  For orthogonal boxes the recip/ucell matrices are diagonal
      // and this reduces to the standard per-axis box subtract.  For non-orthogonal boxes
      // (truncated octahedron etc.) the off-diagonal entries are needed:
      //   df = recip * d  (Cartesian → fractional)
      //   df -= round(df) (wrap into (-0.5, 0.5])
      //   d  = ucell * df (fractional → Cartesian, now the minimum-image displacement)
      if (is_orthog) {
        double dfx = dx * r00;
        double dfy = dy * r11;
        double dfz = dz * r22;
        dfx -= round(dfx);
        dfy -= round(dfy);
        dfz -= round(dfz);
        dx = dfx * u00;
        dy = dfy * u11;
        dz = dfz * u22;
      } else {
        double dfx = dx*r00 + dy*r10 + dz*r20;
        double dfy =          dy*r11 + dz*r21;
        double dfz =                   dz*r22;
        dfx -= round(dfx);
        dfy -= round(dfy);
        dfz -= round(dfz);
        dx = dfx*u00 + dfy*u01 + dfz*u02;
        dy =           dfy*u11 + dfz*u12;
        dz =                     dfz*u22;
      }

      double r2 = dx * dx + dy * dy + dz * dz;

      // Surface marking: if any bath atom is within surface cutoff of this reference atom,
      // mark the ref atom as surface-exposed (-1) so its velocity will also be redirected.
      // For solution regions: surface_cut (3 Å default).  For TM regions: membrane_cut.
      // No race issue: each thread writes only to its own pRecord[orig_r].
      if (r2 < (double)surf_cut2) {
        atomicExch(&pRecord[orig_r], -1);
      }

      // Bath atom group assignment: assign bath atom to the closest reference group within
      // the assignment cutoff (shell_cut for solution, membrane_cut for TM).
      // Pack (float-bits(r2) << 32) | group into 64 bits and atomicMin.  Because positive-float
      // bit patterns sort the same as their float values, atomicMin keeps the smallest distance
      // and its group together — atomically, with no CAS-loop and no (dist, group) tear race.
      // Ties on r2 break toward the lower group number (low bits dominate when high bits match).
      if (r2 < (double)assign_cut2) {
        float r2f = (float)r2;
        unsigned long long packed_new =
            ((unsigned long long)__float_as_uint(r2f) << 32) |
            (unsigned long long)(unsigned int)groupNum;
        atomicMin(&pPacked[orig_s], packed_new);
      }
    }
  }
}

//---------------------------------------------------------------------------------------------
// Double-precision atomicAdd fallback for sm_5x (CC < 6.0).  The native
// atomicAdd(double*, double) is only provided for compute capability >= 6.0, so we supply a
// CAS-loop equivalent on older arches.  Matches the existing convention used in
// kCalculateNEBForces.cu and kNorm.h.
//---------------------------------------------------------------------------------------------
#if (__CUDA_ARCH__ < 600) || (__CUDACC_VER_MAJOR__ < 8)
__device__ void kKAMD_voidatomicAdd(double* ptr, double val) {
  unsigned long long int* ulli_ptr = (unsigned long long int*)ptr;
  unsigned long long int old = *ulli_ptr;
  unsigned long long int assumed;
  do {
    assumed = old;
    old = atomicCAS(ulli_ptr, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (old != assumed);
}
#endif

//---------------------------------------------------------------------------------------------
// kKAMD_ComputeMomenta_kernel
//
// Purpose:
//   Compute the translational momentum P = sum(m_i * v_i) for each reference group entirely
//   on the GPU.  The Fortran side monitors these momenta to detect oscillation completion
//   (2 zero-crossings per group), which triggers the next pseudo-BD reset.
//
//   This kernel replaces downloading the entire velocity array (~706 KB for a 29K-atom system)
//   and computing momenta on the host.  The current host monitor path downloads the first
//   18 doubles (144 bytes), a ~5000x reduction in data transfer.
//
// NL usage:
//   Uses pImageAtomLookup to map each reference atom's original index to its current sorted
//   position, then reads mass and velocity from the image arrays (pImageMass, pImageVelX/Y/Z)
//   in sorted order.
//
// Thread mapping:
//   One thread per reference-group atom.  Each thread contributes its atom's momentum to the
//   group total via atomicAdd on doubles (native on CC >= 6.0; CAS-loop fallback below).
//
// Arguments:
//   pMomenta       - output: per-group momentum laid out as [group*3 + xyz]
//   totalRefAtoms  - total number of reference atoms across all active groups
//   pRefAtomInfo   - per-ref-atom: int2{original_0based_index, group_1toN}
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(256, 1)
kKAMD_ComputeMomenta_kernel(
    double*        pMomenta,
    int            totalRefAtoms,
    int2*          pRefAtomInfo)
{
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= (unsigned int)totalRefAtoms) return;

  int2 refInfo = pRefAtomInfo[tid];
  int orig     = refInfo.x;           // original atom index (0-based)
  int group    = refInfo.y;           // group number (1..ngroups)

  // Map original index -> sorted position, then read mass and velocity
  unsigned int sorted = cSim.pImageAtomLookup[orig];
  double m   = cSim.pImageMass[sorted];
  double vx  = cSim.pImageVelX[sorted];
  double vy  = cSim.pImageVelY[sorted];
  double vz  = cSim.pImageVelZ[sorted];

  // Accumulate momentum for this group (atomic to handle multiple threads per group)
  int base = (group - 1) * 3;
#if (__CUDA_ARCH__ < 600) || (__CUDACC_VER_MAJOR__ < 8)
  kKAMD_voidatomicAdd(&pMomenta[base + 0], m * vx);
  kKAMD_voidatomicAdd(&pMomenta[base + 1], m * vy);
  kKAMD_voidatomicAdd(&pMomenta[base + 2], m * vz);
#else
  atomicAdd(&pMomenta[base + 0], m * vx);
  atomicAdd(&pMomenta[base + 1], m * vy);
  atomicAdd(&pMomenta[base + 2], m * vz);
#endif
}

//---------------------------------------------------------------------------------------------
// kKAMD_ClearMomenta_kernel
//
// Zero the 18-double momenta accumulation buffer before kKAMD_ComputeMomenta_kernel runs.
// Launched with a single warp (1 block, 32 threads); only the first 18 threads write.
//---------------------------------------------------------------------------------------------
__global__ void kKAMD_ClearMomenta_kernel(double* pMomenta)
{
  if (threadIdx.x < 18) pMomenta[threadIdx.x] = 0.0;
}

//---------------------------------------------------------------------------------------------
// kKAMD_ClearBuffers_kernel
//
// Initialize the per-atom result buffers before first-shell detection:
//   pRecord[i] = 0           (no group assigned, not flagged)
//   pPacked[i] = (bits(kamd_cut2) << 32) | 0
//                            (squared cutoff ceiling in the high bits, group=0 in the low bits;
//                             any subsequent atomicMin with a smaller distance will replace it)
//
// Uses a grid-stride loop to cover all atoms regardless of grid size.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(256, 1)
kKAMD_ClearBuffers_kernel(int* pRecord, unsigned long long* pPacked,
                          float kamd_cut2, int natoms)
{
  unsigned long long initPacked = (unsigned long long)__float_as_uint(kamd_cut2) << 32;
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  while (pos < (unsigned int)natoms) {
    pRecord[pos] = 0;
    pPacked[pos] = initPacked;
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kKAMD_ExtractGroup_kernel
//
// Post-pass for kKAMD_FindFirstShell_kernel: copy the group number from the low 32 bits of
// pPacked[i] into pRecord[i], so the redirect kernel can use the existing pRecord[orig] API
// (0 = unaffected, 1..ngroups = solvent group, -1 = surface ref).
//
// Surface-mark sentinel (-1) is preserved: we only overwrite pRecord[i] when it's >= 0.
// (A ref atom marked surface should keep its -1 even if its index also happened to receive
// a group assignment — though for non-overlapping ref/solvent ranges that situation cannot
// arise.)
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(256, 1)
kKAMD_ExtractGroup_kernel(int* pRecord, unsigned long long* pPacked, int natoms)
{
  unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
  while (pos < (unsigned int)natoms) {
    int existing = pRecord[pos];
    if (existing >= 0) {
      int group = (int)(pPacked[pos] & 0xFFFFFFFFu);
      if (group > 0) pRecord[pos] = group;
    }
    pos += blockDim.x * gridDim.x;
  }
}

//---------------------------------------------------------------------------------------------
// kKAMD_RedirectVelocities_kernel
//
// Purpose:
//   Redirect velocities of flagged atoms (first-shell solvent + surface reference atoms)
//   along a random unit vector for their group, while preserving each atom's speed.
//   This creates coherent first-shell motion that accelerates conformational sampling.
//   Operates entirely on GPU velocity buffers — no host round-trip needed.
//
// Leapfrog half-step correction:
//   PMEMD uses the leapfrog integrator, where velocities are stored at half-step (t - dt/2).
//   To redirect the "physical" velocity, we must:
//     1. Advance to full step:    v_full = v_half + F * (dt/2) / m
//     2. Redirect at full step:   v_new  = direction * |v_full|    (preserves speed = KE = T)
//     3. Step back to half step:  v_half_new = v_new - F * (dt/2) / m
//   This ensures the redirected velocity is consistent with the leapfrog time integration.
//
// Force handling:
//   Forces are stored as PMEAccumulator (long long int, fixed-point) and must be descaled
//   by ONEOVERFORCESCALE to get real units.
//   - NTP (constant pressure, non-MPI): bonded forces in pForceXAccumulator, nonbond forces
//     in pNBForceXAccumulator — must sum both.  Controlled by splitForces flag.
//   - NVT / NVE: bonded and nonbond share the same accumulator, so pForceXAccumulator alone
//     contains the total force.  splitForces = 0.
//
// NL usage:
//   This kernel iterates over ALL atoms in sorted order.  For each atom, it uses
//   pImageAtom[sorted_pos] to recover the original atom index, then checks pRecord[orig]
//   to determine whether the atom was flagged by the first-shell kernel.  For surface
//   reference atoms (record == -1), it reads the group assignment from pRefGroupMap[orig].
//   Velocity and force arrays are accessed directly in sorted (image) order.
//
// Thread mapping:
//   One thread per atom (total atoms).  Most threads early-exit (record == 0).
//   Only flagged atoms (~5-10% of the system) do real work.
//
// Arguments:
//   pRecord       - per-atom flag from FindFirstShell (0=skip, 1..ngroups=solvent group, -1=surface)
//   pRandomDir    - ngroups normalized direction vectors: [group*3 + xyz]
//   pRefGroupMap  - per-atom group assignment for reference atoms (0 if not a ref atom)
//   half_dtx      - half of the scaled timestep: 0.5 * dt * sqrt(418.4)
//   natoms        - total number of atoms in the system
//   splitForces   - 1 if NTP (bonded + nonbond in separate accumulators), 0 if combined
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(256, 1)
kKAMD_RedirectVelocities_kernel(
    int*           pRecord,
    double*        pRandomDir,
    int*           pRefGroupMap,
    double         half_dtx,
    int            natoms,
    int            splitForces,
    int            shellMode)
{
  unsigned int sorted_pos = blockIdx.x * blockDim.x + threadIdx.x;
  if (sorted_pos >= (unsigned int)natoms) return;

  // Map sorted position -> original atom index via pImageAtom
  unsigned int orig = cSim.pImageAtom[sorted_pos];
  // Molecular shell mode handles all membrane atoms in a separate one-thread-per-lipid pass.
  if (shellMode && (int)orig >= cSim.KAMDMembraneStart &&
      (int)orig <= cSim.KAMDMembraneEnd) return;
  int record = pRecord[orig];

  if (record == 0) return;  // Not flagged by first-shell detection, skip

  // Determine which group's direction vector to use
  int groupIdx;  // 0-based index into pRandomDir (0..ngroups-1)
  if (record > 0) {
    // Solvent atom assigned to reference group (record = 1..ngroups)
    groupIdx = record - 1;
  } else {
    // Surface reference-group atom (record == -1): look up its group from the static map
    int grp = pRefGroupMap[orig];
    if (grp <= 0) return;  // safety: not actually a reference atom
    groupIdx = grp - 1;
  }

  // Read the pre-normalized random direction for this group (uploaded from host, 144 bytes)
  double dirX = pRandomDir[groupIdx * 3 + 0];
  double dirY = pRandomDir[groupIdx * 3 + 1];
  double dirZ = pRandomDir[groupIdx * 3 + 2];

  // Read current velocity (half-step, stored in image/sorted order)
  double velX = cSim.pImageVelX[sorted_pos];
  double velY = cSim.pImageVelY[sorted_pos];
  double velZ = cSim.pImageVelZ[sorted_pos];

  // Read force from fixed-point accumulators and descale to real units.
  // pForceXAccumulator stores bonded forces (always present).
  // pNBForceXAccumulator stores nonbond forces (separate buffer in NTP, same buffer in NVT).
  PMEAccumulator fx = cSim.pForceXAccumulator[sorted_pos];
  PMEAccumulator fy = cSim.pForceYAccumulator[sorted_pos];
  PMEAccumulator fz = cSim.pForceZAccumulator[sorted_pos];
  if (splitForces) {
    fx += cSim.pNBForceXAccumulator[sorted_pos];
    fy += cSim.pNBForceYAccumulator[sorted_pos];
    fz += cSim.pNBForceZAccumulator[sorted_pos];
  }
  double forceX = (double)fx * (double)ONEOVERFORCESCALE;
  double forceY = (double)fy * (double)ONEOVERFORCESCALE;
  double forceZ = (double)fz * (double)ONEOVERFORCESCALE;

  // wfac = (1/m) * (dt/2) — the velocity increment per unit force over a half-step
  double invMass = cSim.pImageInvMass[sorted_pos];
  double wfac = invMass * half_dtx;

  // Step 1: Advance velocity from half-step to full-step
  double vfx = velX + forceX * wfac;
  double vfy = velY + forceY * wfac;
  double vfz = velZ + forceZ * wfac;

  // Region-type-aware redirect.
  //   Solution (type 1): full 3D redirect along ĝ, preserving total speed.
  //   Transmembrane (type 2): preserve velocity component along the membrane
  //   axis; randomize only the perpendicular plane to a 2D unit direction
  //   (which the host has already projected and renormalized into pRandomDir
  //   for TM groups), preserving the perpendicular speed.
  // ASSUMPTION (pending Haixin confirmation): TM region SURFACE ref atoms get
  // the same membrane-axis-preserving rule as TM first-shell membrane atoms.
  // This mirrors the solution-region symmetry where surface ref atoms get the
  // same redirect as first-shell solvent.  If Haixin says otherwise, this
  // single branch is the only place that needs to change.
  int  regionType = cSim.pKAMDRegionType[groupIdx];
  double vnx, vny, vnz;
  if (regionType == 2) {  // KAMD_RT_TRANSMEMBRANE
    int axis = cSim.KAMDMembraneAxis;  // 1=x, 2=y, 3=z
    double v_par;
    if      (axis == 1) v_par = vfx;
    else if (axis == 2) v_par = vfy;
    else                v_par = vfz;
    double speed_full2 = vfx * vfx + vfy * vfy + vfz * vfz;
    double speed_perp2 = speed_full2 - v_par * v_par;
    if (speed_perp2 < 0.0) speed_perp2 = 0.0;  // guard against FP drift
    double speed_perp  = sqrt(speed_perp2);
    // pRandomDir for TM groups has the parallel component already zeroed and
    // is renormalized to unit length on the perpendicular plane.
    vnx = dirX * speed_perp;
    vny = dirY * speed_perp;
    vnz = dirZ * speed_perp;
    if      (axis == 1) vnx += v_par;
    else if (axis == 2) vny += v_par;
    else                vnz += v_par;
  } else {                // KAMD_RT_SOLUTION
    double speed = sqrt(vfx * vfx + vfy * vfy + vfz * vfz);
    vnx = dirX * speed;
    vny = dirY * speed;
    vnz = dirZ * speed;
  }

  // Step back to half-step for leapfrog consistency.
  cSim.pImageVelX[sorted_pos] = vnx - forceX * wfac;
  cSim.pImageVelY[sorted_pos] = vny - forceY * wfac;
  cSim.pImageVelZ[sorted_pos] = vnz - forceZ * wfac;
}

// One thread owns one complete lipid residue, so no atom can be written twice.
// The first-shell atom assignments are converted into a mass-weighted molecular
// owner. Only the molecule's lateral COM velocity is turned; internal velocities
// and the component along the membrane normal are exactly preserved.
__global__ void
kKAMD_RedirectMembraneMolecules_kernel(int2* molRange, int* previousOwner,
    int nMol, int* record, double* randomDir, int nGroups,
    double half_dtx, int splitForces, double maxAngle, double* impulse)
{
  int m = blockIdx.x * blockDim.x + threadIdx.x;
  if (m >= nMol) return;
  int2 range = molRange[m];
  double score[32];
  for (int g = 0; g < 32; g++) score[g] = 0.0;
  double massTotal = 0.0;
  double momentum[3] = {0.0, 0.0, 0.0};

  for (int orig = range.x; orig <= range.y; orig++) {
    unsigned int s = cSim.pImageAtomLookup[orig];
    double mass = cSim.pImageMass[s];
    int group = record[orig];
    if (group > 0 && group <= nGroups && cSim.pKAMDRegionType[group-1] == 2)
      score[group-1] += mass;
    massTotal += mass;
    PMEAccumulator fx = cSim.pForceXAccumulator[s];
    PMEAccumulator fy = cSim.pForceYAccumulator[s];
    PMEAccumulator fz = cSim.pForceZAccumulator[s];
    if (splitForces) {
      fx += cSim.pNBForceXAccumulator[s];
      fy += cSim.pNBForceYAccumulator[s];
      fz += cSim.pNBForceZAccumulator[s];
    }
    double wfac = cSim.pImageInvMass[s] * half_dtx * (double)ONEOVERFORCESCALE;
    momentum[0] += mass * (cSim.pImageVelX[s] + (double)fx * wfac);
    momentum[1] += mass * (cSim.pImageVelY[s] + (double)fy * wfac);
    momentum[2] += mass * (cSim.pImageVelZ[s] + (double)fz * wfac);
  }
  // Include every lipid in the denominator, even when it does not contact a
  // reference group. This makes the correction a uniform membrane COM shift.
  atomicAdd(&impulse[3], massTotal);

  int owner = 0;
  double best = 0.0;
  for (int g = 0; g < nGroups; g++) {
    if (score[g] > best) { best = score[g]; owner = g + 1; }
  }
  if (best <= 0.0 || massTotal <= 0.0) { previousOwner[m] = 0; return; }
  int old = previousOwner[m];
  // Keep an existing owner while its contact mass remains competitive.
  if (old > 0 && old <= nGroups && score[old-1] >= 0.8 * best) owner = old;
  previousOwner[m] = owner;

  int axis = cSim.KAMDMembraneAxis;
  int a = axis == 1 ? 1 : 0;
  int b = axis == 3 ? 1 : 2;
  double vx = momentum[a] / massTotal;
  double vy = momentum[b] / massTotal;
  double speed = sqrt(vx*vx + vy*vy);
  if (speed < 1.0e-12) return;
  double dx = randomDir[3*(owner-1) + a];
  double dy = randomDir[3*(owner-1) + b];
  double angle = atan2(vx*dy - vy*dx, vx*dx + vy*dy);
  // Weak contact gives a smaller turn; a complete shell contact gets the full cap.
  double cap = maxAngle * sqrt(score[owner-1] / massTotal);
  if (angle > cap) angle = cap;
  if (angle < -cap) angle = -cap;
  double cs = cos(angle), sn = sin(angle);
  double deltaA = vx*cs - vy*sn - vx;
  double deltaB = vx*sn + vy*cs - vy;
  atomicAdd(&impulse[a], massTotal * deltaA);
  atomicAdd(&impulse[b], massTotal * deltaB);
  atomicAdd(&impulse[4], 1.0);
  for (int orig = range.x; orig <= range.y; orig++) {
    unsigned int s = cSim.pImageAtomLookup[orig];
    if (a == 0) cSim.pImageVelX[s] += deltaA;
    else        cSim.pImageVelY[s] += deltaA;
    if (b == 1) cSim.pImageVelY[s] += deltaB;
    else        cSim.pImageVelZ[s] += deltaB;
  }
}

// Remove only the net in-plane momentum that molecular shell redirection added.
// All lipids receive the same small velocity offset, so the local, coherent
// first-shell motion remains available to make space for the moving protein.
__global__ void
kKAMD_BalanceMembraneImpulse_kernel(int first, int nAtoms, int axis,
                                    const double* impulse)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= nAtoms || impulse[3] <= 0.0) return;
  int orig = first + index;
  unsigned int s = cSim.pImageAtomLookup[orig];
  double invMass = 1.0 / impulse[3];
  if (axis != 1) cSim.pImageVelX[s] -= impulse[0] * invMass;
  if (axis != 2) cSim.pImageVelY[s] -= impulse[1] * invMass;
  if (axis != 3) cSim.pImageVelZ[s] -= impulse[2] * invMass;
}

//=============================================================================================
// Host-side launchers
//
// These functions are called from the extern "C" wrappers in gpu.cpp, which in turn are
// called from the Fortran MD loop in runmd.F90.  They configure grid dimensions and launch
// the device kernels above.
//=============================================================================================

//---------------------------------------------------------------------------------------------
// kKAMD_FindFirstShell: clear result buffers, then launch the first-shell detection kernel,
// then extract the per-atom group from the packed buffer into pRecord.
//
// Sequence:
//   1. Clear pRecord (to 0) and pPacked (to bits(kamd_cut2) << 32) for all atoms.
//   2. Launch the search kernel with one thread per reference atom.  Surface marks land in
//      pRecord directly; group assignments race-safely into pPacked via atomicMin.
//   3. Extract the group number from pPacked into pRecord (preserving -1 surface marks).
//
// The isNTP flag is set for constant-pressure simulations (ntp > 0) EXCEPT Monte Carlo
// barostat (barostat == 2), which does not update ucell continuously.
//---------------------------------------------------------------------------------------------
extern "C" void kKAMD_FindFirstShell(gpuContext gpu, int totalRefAtoms)
{
  // Initialize pPacked with a huge value in the high 32 bits.  Any subsequent
  // atomicMin from a real (dist²,group) candidate will replace it; this saves
  // having to pick an "outer" cutoff that covers both region types.
  float initCut2 = 1.0e30f;
  int clearBlocks = (gpu->sim.atoms + 255) / 256;
  kKAMD_ClearBuffers_kernel<<<clearBlocks, 256>>>(
      gpu->sim.pKAMDRecord, gpu->sim.pKAMDPacked, initCut2, gpu->sim.atoms);
  LAUNCHERROR("kKAMD_ClearBuffers");

  if (totalRefAtoms > 0) {
    int isNTP = 0;
    if ((gpu->sim.ntp > 0) && (gpu->sim.barostat != 2)) {
      isNTP = 1;
    }
    int searchBlocks = (totalRefAtoms + 255) / 256;
    kKAMD_FindFirstShell_kernel<<<searchBlocks, 256>>>(
        gpu->sim.pKAMDRecord, gpu->sim.pKAMDPacked,
        totalRefAtoms,
        (int2*)gpu->pbKAMDRefAtomInfo->_pDevData,
        isNTP);
    LAUNCHERROR("kKAMD_FindFirstShell");

    kKAMD_ExtractGroup_kernel<<<clearBlocks, 256>>>(
        gpu->sim.pKAMDRecord, gpu->sim.pKAMDPacked, gpu->sim.atoms);
    LAUNCHERROR("kKAMD_ExtractGroup");
  }
}

//---------------------------------------------------------------------------------------------
// kKAMD_RedirectVelocities: launch the velocity redirection kernel for all flagged atoms.
//
// The splitForces flag mirrors the logic in kForcesUpdate.cu: in single-GPU NTP simulations,
// bonded and nonbond forces are accumulated in separate buffers (pForceXAccumulator and
// pNBForceXAccumulator).  In MPI or NVT runs, they share the same buffer.
//---------------------------------------------------------------------------------------------
extern "C" void kKAMD_RedirectVelocities(gpuContext gpu,
                                           double* pRandomDir,
                                           int* pRefGroupMap,
                                           double half_dtx)
{
  int splitForces = 0;
#ifndef MPI
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat != 2)) {
    splitForces = 1;
  }
#endif

  int blocks = (gpu->sim.atoms + 255) / 256;
  kKAMD_RedirectVelocities_kernel<<<blocks, 256>>>(
      gpu->sim.pKAMDRecord,
      pRandomDir,
      pRefGroupMap,
      half_dtx,
      gpu->sim.atoms,
      splitForces,
      gpu->kamdShellMode);
  LAUNCHERROR("kKAMD_RedirectVelocities");
}

extern "C" void kKAMD_RedirectMembraneMolecules(gpuContext gpu,
                                                  double* pRandomDir,
                                                  double half_dtx)
{
  if (gpu->kamdShellMode == 0 || gpu->nKAMDMembraneMol == 0) return;
  cudaMemsetAsync(gpu->pbKAMDMembraneImpulse->_pDevData, 0, 5 * sizeof(double));
  LAUNCHERROR("kKAMD_MembraneImpulse_zero");
  int splitForces = 0;
#ifndef MPI
  if ((gpu->sim.ntp > 0) && (gpu->sim.barostat != 2)) splitForces = 1;
#endif
  int blocks = (gpu->nKAMDMembraneMol + 127) / 128;
  kKAMD_RedirectMembraneMolecules_kernel<<<blocks, 128>>>(
      gpu->pbKAMDMembraneMol->_pDevData,
      gpu->pbKAMDMembraneOwner->_pDevData,
      gpu->nKAMDMembraneMol,
      gpu->sim.pKAMDRecord,
      pRandomDir,
      gpu->pbKAMDRegionType->_length,
      half_dtx,
      splitForces,
      gpu->kamdShellMaxAngle,
      gpu->pbKAMDMembraneImpulse->_pDevData);
  LAUNCHERROR("kKAMD_RedirectMembraneMolecules");
  int nAtoms = gpu->sim.KAMDMembraneEnd - gpu->sim.KAMDMembraneStart + 1;
  int balanceBlocks = (nAtoms + 255) / 256;
  kKAMD_BalanceMembraneImpulse_kernel<<<balanceBlocks, 256>>>(
      gpu->sim.KAMDMembraneStart, nAtoms, gpu->sim.KAMDMembraneAxis,
      gpu->pbKAMDMembraneImpulse->_pDevData);
  LAUNCHERROR("kKAMD_BalanceMembraneImpulse");
}

//---------------------------------------------------------------------------------------------
// kKAMD_ComputeMomenta: zero the accumulation buffer, then compute per-group momenta.
//
// Two-kernel sequence:
//   1. kKAMD_ClearMomenta_kernel — zero 18 doubles (single warp).
//   2. kKAMD_ComputeMomenta_kernel — one thread per ref atom, atomicAdd contributions.
// The caller (gpu_kamd_compute_momenta_ in gpu.cpp) issues cudaDeviceSynchronize and then
// downloads the 18 doubles to the host for the Fortran momentum-monitoring logic.
//---------------------------------------------------------------------------------------------
extern "C" void kKAMD_ComputeMomenta(gpuContext gpu, int totalRefAtoms)
{
  // Zero the entire ngroups*3 doubles momenta buffer (was hardcoded to 18 doubles =
  // first 6 groups, a leftover from the pre-namelist 6-group limit).
  cudaMemsetAsync(gpu->pbKAMDMomenta->_pDevData, 0,
                  gpu->pbKAMDMomenta->_length * sizeof(double));
  LAUNCHERROR("kKAMD_ClearMomenta");

  if (totalRefAtoms > 0) {
    int blocks = (totalRefAtoms + 255) / 256;
    kKAMD_ComputeMomenta_kernel<<<blocks, 256>>>(
        gpu->pbKAMDMomenta->_pDevData,
        totalRefAtoms,
        (int2*)gpu->pbKAMDRefAtomInfo->_pDevData);
    LAUNCHERROR("kKAMD_ComputeMomenta");
  }
}

//---------------------------------------------------------------------------------------------
// kKAMD_GatherIons_kernel
//   One thread per ion: map original ion index -> sorted position, read xyz, write to a
//   compact output buffer in ion order.  No atomics (each thread owns its slot).
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(256, 1)
kKAMD_GatherIons_kernel(double* pIonCrd, int nIons, int* pIonAtoms)
{
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= (unsigned int)nIons) return;

  unsigned int sorted = cSim.pImageAtomLookup[pIonAtoms[tid]];
  pIonCrd[tid * 3 + 0] = cSim.pImageX[sorted];
  pIonCrd[tid * 3 + 1] = cSim.pImageY[sorted];
  pIonCrd[tid * 3 + 2] = cSim.pImageZ[sorted];
}

//---------------------------------------------------------------------------------------------
// kKAMD_MemCOM_kernel
//   Grid-stride over the contiguous membrane range [mStart0, mEnd0] (0-based, inclusive);
//   mass-weighted coordinate accumulation into pCOMacc = {mx, my, mz, m} via atomicAdd.
//   Caller zeros pCOMacc before launch; host divides {mx,my,mz} by m to obtain the COM.
//---------------------------------------------------------------------------------------------
__global__ void
__LAUNCH_BOUNDS__(256, 1)
kKAMD_MemCOM_kernel(double* pCOMacc, int mStart0, int mEnd0)
{
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int stride = blockDim.x * gridDim.x;
  for (int i = mStart0 + (int)tid; i <= mEnd0; i += (int)stride) {
    double m = cSim.pImageMass[i];
#if (__CUDA_ARCH__ < 600) || (__CUDACC_VER_MAJOR__ < 8)
    kKAMD_voidatomicAdd(&pCOMacc[0], m * cSim.pImageX[i]);
    kKAMD_voidatomicAdd(&pCOMacc[1], m * cSim.pImageY[i]);
    kKAMD_voidatomicAdd(&pCOMacc[2], m * cSim.pImageZ[i]);
    kKAMD_voidatomicAdd(&pCOMacc[3], m);
#else
    atomicAdd(&pCOMacc[0], m * cSim.pImageX[i]);
    atomicAdd(&pCOMacc[1], m * cSim.pImageY[i]);
    atomicAdd(&pCOMacc[2], m * cSim.pImageZ[i]);
    atomicAdd(&pCOMacc[3], m);
#endif
  }
}

//---------------------------------------------------------------------------------------------
// kKAMD_GatherIons / kKAMD_MemCOM launchers.
//---------------------------------------------------------------------------------------------
extern "C" void kKAMD_GatherIons(gpuContext gpu, int nIons)
{
  if (nIons <= 0) return;
  int blocks = (nIons + 255) / 256;
  kKAMD_GatherIons_kernel<<<blocks, 256>>>(
      gpu->pbKAMDIonCrd->_pDevData, nIons, gpu->pbKAMDIons->_pDevData);
  LAUNCHERROR("kKAMD_GatherIons");
}

extern "C" void kKAMD_MemCOM(gpuContext gpu, int mStart0, int mEnd0)
{
  int nAtoms = mEnd0 - mStart0 + 1;
  if (nAtoms <= 0) return;
  int blocks = (nAtoms + 255) / 256;
  kKAMD_MemCOM_kernel<<<blocks, 256>>>(
      gpu->pbKAMDIonCOMacc->_pDevData, mStart0, mEnd0);
  LAUNCHERROR("kKAMD_MemCOM");
}
