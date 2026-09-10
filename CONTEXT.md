# Repository context

Last checked: 2026-09-10. Read this file when resuming work on PfQMC.
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
  Its source derivation is section 3.2 of the external notes under
  `/home/han/notes/sign_free_classification/`; those notes are not bundled here.

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
