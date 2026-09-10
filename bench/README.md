# PFQMC performance experiments

This worktree accelerates the existing signed-Pfaffian Monte Carlo implementation.
The main improvement applies the interaction layer's disjoint 2x2 blocks directly
instead of using dense matrix products. It applies to both sweep directions,
segment construction, and interaction multiplication in raw-sign initialization.
The UDT permutation is also formed by directly scattering triangular columns.
The RNG sequence, proposal order, acceptance phase, kinetic operator, QR pivoting,
stabilization interval, and signed Pfaffian evaluation are retained.

## Reproduce

Run in WSL Ubuntu as `han`, using the oneAPI environment in [CONTEXT.md](../CONTEXT.md).
All commands below are from the worktree root. Eigen and CMake are reused from
the original checkout; builds and outputs are local to this worktree.

```bash
source /opt/intel/oneapi/setvars.sh
export PATH="/home/han/PfQMC/build/python-tools/cmake/data/bin:$PATH"
make -C inc/pfapack/fortran mFC=ifx
make -C inc/pfapack/c_interface
cmake -S . -B build/perf -DCMAKE_CXX_COMPILER=mpiicpx \
  -DEIGEN3_INCLUDE_DIR=/home/han/PfQMC/build/eigen \
  -DCMAKE_PREFIX_PATH=/home/han/miniconda3
cmake --build build/perf -j 2
ctest --test-dir build/perf --output-on-failure
cp build/perf/pfqmc_bench build/perf/pfqmc_bench_blocks
/home/han/miniconda3/bin/python bench/build_reference.py
/home/han/miniconda3/bin/python bench/benchmark.py \
  --output build/perf/reproduction.csv --versions baseline blocks \
  --suite large --repeats 3
/home/han/miniconda3/bin/python bench/summarize.py build/perf/reproduction.csv
```

`build_reference.py` copies the current build configuration, dependencies, and
benchmark harness under ignored `build/`, restores `inc/qr_udt.h` and
`inc/spinless_tV.h` from `3832a1f8e9efb8cba24b52a58ddc843bfc4beb7f`, and builds
`pfqmc_bench_baseline`. These are the only production files changed by the final
performance implementation. It never edits the active checkout's source files.
It requires PFAPACK's archives to have already been built as above. The inherited
WSL build/documentation repairs are recorded separately in commit `ec09e0c`.

For a short pilot comparison use `--suite pilot --repeats 1`.
For a single run:

```bash
build/perf/pfqmc_bench majorana 12 12 8 1 42 1
# MODEL L Ltau timed_sweep_pairs MKL_threads seed validate
```

Models are `majorana`, `honeycomb`, `square`, `chain0`, and `chain1`. For 2D models,
L means Lx=Ly; for chains it means sites. Both honeycomb models have 2*L*L physical
sites, but the full honeycomb matrix has twice the single-Majorana dimension.
Square/chain matrices have twice as many Majorana modes as physical sites.
Use even L for the square and chain benchmarks.

## Measurement method

All timings use `std::chrono::steady_clock` and one process. MKL dynamic threading
is disabled; both binaries use identical explicit thread counts and compiler flags.
The physical parameters are dt=0.1, V=0.7, stabilization interval=10; chains also
use delta=0.3, mu=0.2 and periodic boundaries. Model construction and initial raw
signed Pfaffian calculation are measured as `init_s`. One right+left sweep pair
warms the code; then the requested complete pairs are measured as `sweeps_s`.
`total_s` includes initialization, warmup and timed sweeps. Post-timing validation
and output serialization are excluded. This measures the actual sampler's sweep
throughput, not a stand-alone matrix microbenchmark or full CLI file-I/O cost.

Baseline and optimized processes run sequentially, alternating their order on
successive repeats. No compilation or other benchmark runs concurrently. The
runner compares final auxiliary-field checksums, Green-matrix checksums, and
complex signed phases across versions and fails on disagreement. Each binary
also checks the sign against `getSignRaw()` and rebuilds Green's matrix. These are
numerical consistency checks, not a proof of sign freedom or an autocorrelation /
effective-independent-sample study. Identical full auxiliary configurations and
Green matrices are additionally checked after each sweep by the regression tests.

## Iterations

1. **UDT column scatter:** form D^-1 R P^T by scattering the upper-triangular
   columns using LAPACK's pivot indices. This replaces an O(N^3) triangular-times-
   permutation multiply with O(N^2) writes. Whole-sweep differences in the first
   three-repeat pilot were only about 0-2%, sometimes within run-to-run noise.
   This is a small structural simplification, not the source of the headline gain.
2. **Interaction blocks:** detect fixed nonzero pairs once, then apply each 2x2
   block directly in O(N^2) instead of O(N^3). Reads the updated matrix entries
   after accepted flips. Supports in-place multiplication and identity rows at
   open boundaries. Overlapping supports or unpaired nonidentity diagonals fall
   back to the original dense product. This is the main measured improvement.
3. **Fused rank-two update (rejected):** combine the two unconjugated BLAS GERU
   traversals in an Eigen expression. Numerical tests passed, but the sweep was
   51-58% slower than the block version in the pilot. The final production code
   retains MKL GERU. [The experimental patch](experiments/fused_rank2.patch) is kept
   only for reproduction; it is not applied. From the final source, `git apply`
   of that patch recreates the candidate, then build/copy the benchmark as
   `build/perf/pfqmc_bench_fused` and compare with `--versions blocks fused`.
   [The extra candidate-only regression](experiments/fused_rank2_test.cpp) can be
   appended to `test/performanceTest.cpp` when reproducing its 19-test run.

Pilot raw measurements: [stage 1](results/stage1.csv), [stage 2](results/stage2.csv),
[stage 3](results/stage3.csv). Final measurements and machine details are recorded
in [RESULTS.md](RESULTS.md).

## Regression coverage and scope

The added tests compare block operations against dense matrix expressions before
and after accepted flips, including aliasing, both propagation directions,
stabilized products, honeycomb, single Majorana, square, both chain HS channels,
open boundaries, overlapping odd-periodic supports, and V=0. Seeded trajectories
are compared with a dense interaction oracle, including incomplete stabilization
segments. A synthetic negative-ratio case explicitly checks that the signed phase
is retained. A scaled-column matrix checks UDT reconstruction and Q unitarity.
The original tests and the four-mode exact Fock/signed-Pfaffian validator also run.

This does not make previously unsupported odd lattices or extreme imaginary-time
regimes valid. The dense fallback preserves existing behavior for overlapping
layers. Results describe the tested oneAPI/WSL laptop and parameter ranges;
multinode MPI and GPU execution have not been benchmarked. Dense kinetic products,
pivoted QR, Green reconstruction, raw Pfaffians, and immediate accepted updates
remain; the complete algorithm still has cubic work as lattice size grows.
