# Repository context

Last checked: 2026-09-11. Read this file when resuming work on PfQMC.
This describes the verified local machine, not a portable installation guarantee.

## Working rules and project map

- All repository file and directory names must use English/ASCII names. Document
  prose may be Chinese or English. Keep links, CMake paths, commands, and ignore
  rules consistent when renaming anything.
- Start with `git status --short` and `git log -3 --oneline`; preserve existing
  user changes. Update this context when the environment or workflow changes.
- Run builds, Git, and executables inside WSL Ubuntu as user `han`.
- `main.cpp`: CLI and original square, honeycomb, Majorana honeycomb, and Kitaev
  chain simulations. `inc/`: model/operator headers. `src/`: PfQMC and skew-matrix
  utilities. `test/`: original GoogleTest regressions.
- `inc/pfapack/`: bundled PFAPACK Fortran routines and C interface.
- `sign_free_verify/four_mode.cpp`: independent four-complex-mode sign validation;
  CMake executable target is **`sign_free_check`**, CTest name is
  **`SignFree.FourMode`**. Source directory and executable intentionally differ.
- [Model and reproduction](sign_free_verify/README.md) and
  [historical results](sign_free_verify/REPORT.md) document the four-mode model.
  Its source derivation is section 3.2 of the external notes now under
  `/home/han/notes/sign_free/sign_free_classification/`; those notes are not bundled here.
- `inc/catalog/`, `src/catalog_*.cpp`, `sign_free_catalog/`: parameterized
  Hamiltonians from the audited sign-free notes, a dense continuous-HS walker
  using existing `PfQMC` sweeps, and independent signed validation. CMake
  library `pfqmc_catalog`, executable `sign_free_models`. See the
  [catalog workflow](sign_free_catalog/README.md), [mathematical models](sign_free_catalog/MODELS.md)
  and [recorded results](sign_free_catalog/RESULTS.md). Static sectors,
  non-Hermitian generators and worldline projectors have separate APIs;
  they must not be treated as ordinary Hermitian Gaussian-HS models.
- Continuing sign-free search: `sign_free_catalog/MODELS.md` section 12,
  `search_registry.json`, `SEARCH_RESULTS.md`, and `SEARCH_QUEUE.md` record
  new constructions, deduplication, signed evidence and unvalidated leads.
  CMake target `sign_free_search` uses `src/catalog_search.cpp`; the old
  `sign_free_models` audit remains a separate historical validation suite.
  Count only newly validated Hamiltonian constructions, not old catalog entries
  or size/parameter variants. Stop expanding after the count exceeds 20, then
  perform the final deduplication/merging and write complete per-family proofs.
  Seven search rounds have 21 distinct validated constructions; the search
  is complete after merging the Ising/XXZ faces; consult the registry
  rather than this snapshot for the current count. Dedicated targets
  `sign_free_projected` and `sign_free_dynamic_ising` validate sector-projected
  and CT-INT configuration weights; they are not production Monte Carlo walkers.
  `sign_free_poisson` validates local exponential CT vertices;
  `sign_free_albert` runs the full Albert Jordan family with actual PfQMC sweeps.
  `sign_free_mobile_spin` and `sign_free_quantum_link` validate dynamic
  composite-pair and one-dimensional correlated-hopping CT paths.
  `sign_free_meron`, `sign_free_ring_exchange`, `sign_free_color_projector`,
  and `sign_free_frustrated_impurity` validate further CT decompositions.
  Further signed CT targets: `sign_free_stabilizer`, `sign_free_jq`,
  `sign_free_bilayer`, `sign_free_zigzag`, `sign_free_tj_chain`, `sign_free_fredkin`.
  `sign_free_pxp` adds full-space kinetic-constraint validation.
  Each round preserves its evidence separately; 84 CTest cases passed.
  `summarize_search_all.py` aggregates seven rounds without mutating the registry
  or historical evidence. MODELS.md has the final 21-family index and proofs.
  Representation correction: the 21-family total is NOT an auxiliary-field PfQMC
  count. Only R001-R005 and R009 (6 families) have confirmed ordinary Hermitian
  Gaussian-HS weights and actual PfQMC sweeps; 15 remain unconfirmed in that
  representation. R006 uses HS plus conserved-sector sums; R008 uses CT Gaussian
  vertices; the others use worldlines/CT/SSE. Legacy sections 8-11 are separately
  qualified (static backgrounds, non-Gaussian unitary mapping, worldlines, and
  formal non-Hermitian Gaussian generators). See MODELS.md's status table and
  each registry entry's auxiliary_field_pfqmc field. Historical evidence stays intact.
  Catalog raw CSVs and per-run logs are local generated artifacts ignored by Git;
  CTest logs, summary JSON, checksum manifests, registry, and small independent
  research fixtures remain eligible for version control. A fresh clone needs
  the raw samples regenerated or restored before running the summary scripts.
  Historical checksum manifests describe the recorded runs, not later edits.

## Verified WSL2 and Intel oneAPI environment

| Component | Local value |
| --- | --- |
| WSL distribution/user | `Ubuntu` / `han` |
| Linux | Ubuntu 24.04.3 LTS, x86_64 |
| Kernel | `6.6.87.2-microsoft-standard-WSL2` |
| Repository in Linux | `/home/han/PfQMC` |
| Repository from Windows | `\\wsl.localhost\Ubuntu\home\han\PfQMC` |
| oneAPI activation | `/opt/intel/oneapi/setvars.sh` |
| C++/C/Fortran | `icpx` / `icx` / `ifx`, compiler 2025.3.2 |
| Compiler directory | `/opt/intel/oneapi/compiler/2025.3/bin` |
| MPI | Intel MPI 2021.17, `/opt/intel/oneapi/mpi/2021.17` |
| MPI C++ wrapper | `/opt/intel/oneapi/mpi/2021.17/bin/mpiicpx` |
| oneMKL | 2025.3.0, `MKLROOT=/opt/intel/oneapi/mkl/2025.3` |
| MKL CMake package | `/opt/intel/oneapi/mkl/2025.3/lib/cmake/mkl` |
| MKL interface/threading | LP64, dynamic libraries, Intel threading |
| Eigen | 3.4.0, `/home/han/PfQMC/build/eigen` (contains `Eigen/`) |
| GoogleTest | 1.11.0, `/home/han/miniconda3/lib/cmake/GTest` |
| CMake / CTest | 4.4.3, `/home/han/PfQMC/build/python-tools/cmake/data/bin` |
| Python for local tool installation | `/home/han/miniconda3/bin/python` |

CMake and Eigen already exist under ignored `build/`; GoogleTest is supplied by
Miniconda. Neither `cmake` nor `ctest` is on the default WSL PATH on this machine.
Source oneAPI and add the CMake directory in **each new shell**, including shells
used only to run tests. Activation supplies compiler/MPI paths and runtime library
search paths. Do not substitute Windows MKL libraries for these Linux libraries.

From PowerShell, enter the Linux environment first:

```powershell
wsl -d Ubuntu -u han
```

Then use Bash:

```bash
cd /home/han/PfQMC
source /opt/intel/oneapi/setvars.sh
export PATH="$PWD/build/python-tools/cmake/data/bin:$PATH"
command -v mpiicpx icx ifx cmake ctest
printf '%s\n' "$MKLROOT"
```

## Build and test

Run from the repository root after the activation above:

```bash
make -C inc/pfapack/fortran mFC=ifx
make -C inc/pfapack/c_interface
cmake -S . -B build/context-check \
  -DCMAKE_CXX_COMPILER=mpiicpx \
  -DEIGEN3_INCLUDE_DIR="$PWD/build/eigen" \
  -DCMAKE_PREFIX_PATH=/home/han/miniconda3
cmake --build build/context-check -j 2
ctest --test-dir build/context-check --output-on-failure
```

`build/context-check` is a separate build tree used for this repair. It avoids
reusing the old `build/CMakeCache.txt`, whose detected compiler was `/usr/bin/c++`
even though the project later selected `mpiicpx`. Set the compiler on the first
configure so detection and compilation agree. The repository's CMake still
selects `mpiicpx` and uses Intel-specific floating-point flags; a generic GCC
build is not the documented workflow. Existing `-ffp-contract` override warnings
have been observed and should not be confused with a failed build.

PFAPACK is built separately because CMake links its two static archives directly;
it does not build them. The C interface uses `icx`, and `mFC=ifx` overrides the
Fortran makefile's legacy `ifort` default. The top-level legacy `Makefile` uses
`mpiicpc` and an AVX-512 target; use the CMake workflow above on this machine.

If moving to a fresh checkout, the local CMake/Eigen directories must first be
restored. These commands download dependencies and are unnecessary when present:

```bash
mkdir -p build
/home/han/miniconda3/bin/python -m pip install --target build/python-tools cmake==4.4.3
git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git build/eigen
```

The commands assume oneAPI and the listed Miniconda GoogleTest package are already
installed. On another machine, adjust the paths. Do not delete all of `build/`
to reset CMake: it also holds the local CMake installation and Eigen source.

## Run the four-mode validator

```bash
build/context-check/sign_free_check --samples 24 --seed 240810311 \
  --csv sign_free_verify/samples_smoke.csv
# Longer reproduction of the historical sample counts:
build/context-check/sign_free_check --samples 10000 --seed 240810311 \
  --csv sign_free_verify/samples_240810311.csv
build/context-check/sign_free_check --samples 1000 --seed 910 \
  --csv sign_free_verify/samples_910.csv
```

`--samples N` gives N configurations in each of four random ensembles plus three
boundary cases: `4*N+3` total. The validator also checks algebra, Gaussian HS
normalization, pairing contributions, and a negative-weight control. It forces
MKL to one thread. All failures return a nonzero exit code.

This is a configuration-weight validator, not a production Monte Carlo walker
for the new model. It compares exact 16-dimensional Fock traces, a 3-dimensional
auxiliary character formula, and signed Pfaffians. The `Z >= 16` bound refers to
the fermion trace before positive HS density/constant factors. Never replace the
signed weight with an absolute value or a positive determinant square root.
Direct products/inverses are only validated over the documented small-system
parameter range; there is no long-imaginary-time stabilization in this validator.

## Original-program smoke runs

After activation and build, from the repository root:

```bash
mkdir -p build/smoke
build/context-check/main --MRhoneycomb 2 2 4 0.1 0.7 1 42 10 build/smoke/
build/context-check/main --chain 4 4 0.1 0.7 0.3 1 42 10 build/smoke/ 0.2 0 1
```

The positional arguments are:

```text
--MRhoneycomb Lx Ly Ltau dt V nthreads seed measurements output_directory
--chain Lx Ltau dt V delta nthreads seed measurements output_directory mu hsScheme boundary
```

Supply every argument, keep the output path short, and include its trailing slash:
the current CLI does not fully validate argument counts and uses fixed-size output
buffers. `--square` and `--honeycomb` run larger hard-coded defaults and may take
several minutes. The square `raw sign ave` output is a placeholder initialized to
1; it is not an independent sign check. See the historical report for details.

These examples run one MPI process. `OFI endpoint open failed` has previously
occurred under a restricted execution sandbox; investigate MPI/network execution
permissions before attributing that startup failure to the numerical algorithm.
Multi-node execution has not been validated by these smoke tests.

## Repair of commit 3832a1f

The latest commit when this task started was `3832a1f` (`sign free verify`).
Its CMake target referred to nonexistent `sign_free_check/four_mode.cpp`, which
caused a fresh configure to fail. The actual source is under `sign_free_verify/`.
The root README links and validator CSV examples also retained an obsolete
non-English directory name. These references now use `sign_free_verify/`.
The CI source-path filter now includes changes to the validator.

The root build example also needed correction: it changed into the PFAPACK
subdirectory without returning to the root, activated oneAPI too late, and used
`&` where sequential configure/build execution was required. It now uses `make -C`
and separate CMake commands with activation first.

The historical report describes earlier local runs; it is not evidence that the
broken committed CMake paths worked. CSV/log/checksum files in `sign_free_verify/`
are ignored, exist locally, and are absent from a fresh clone. Reproduction commands
can overwrite those sample CSVs; preserve them first if comparing historical data.
Record fresh verification separately rather than silently changing historical results.


## Performance worktree and repeatable timing (2026-09-10)

- Independent worktree: `/home/han/PfQMC-perf`, branch `codex/pfqmc-performance`.
  The original checkout's uncommitted environment/documentation repairs were
  copied and recorded separately as `ec09e0c`; original files were not modified.
- This worktree reuses Eigen and CMake from `/home/han/PfQMC/build/` using absolute
  paths; its own generated artifacts are under `/home/han/PfQMC-perf/build/perf`.
  Each build/test/run shell still activates `/opt/intel/oneapi/setvars.sh`.
- [Benchmark workflow](bench/README.md) documents the `pfqmc_bench` target,
  reference-build script, CSV schema and timing scope. The reference uses the
  original kernels from `3832a1f` with the same current benchmark and build flags.
- [Measured results](bench/RESULTS.md) include raw measurements for sites >100,
  retained/rejected iterations, signed-phase checks, and regression results.
  Sweep throughput is measured at fixed physical parameters, RNG seeds, thread
  counts and stabilization schedule; numerical agreement is not a proof.
- WSL reports an Intel Core i9-13900H and 20 logical CPUs. Final benchmarks run
  one process at a time with MKL dynamic threading disabled and explicit 1 or 4
  threads. Do not run builds/tests concurrently with timed comparisons.
- `rg` is unavailable in this WSL installation; use `git grep` or `find`/`grep`.

## Git line endings (2026-09-11)

- WSL Git inherits `core.autocrlf=input` and `core.safecrlf=true` from
  `/home/han/.gitconfig`. Keep repository text files saved with LF line endings.
- Mixed CRLF/LF in `sign_free_catalog/research/dynamic_ising_screen.py`
  blocked `git add` with `CRLF would be replaced by LF`; only its CRLF line
  endings were normalized to LF. Keep the safety check enabled.
- Verified `git add --all` using a temporary index; the real staging area and
  commit history were left unchanged.
