# Large-lattice performance results

Measured on 2026-09-10 in WSL Ubuntu, Intel Core i9-13900H (20 logical CPUs),
oneAPI C++ 2025.3.2 / oneMKL 2025.3.0, Eigen 3.4.0. Both versions use the same
`mpiicpx` compiler and `-fp-model=precise -fimf-arch-consistency=true -no-fma -O3`.
One process runs at a time; physical parameters and seeds are held fixed.

The baseline contains kernels from `3832a1f8e9efb8cba24b52a58ddc843bfc4beb7f` plus
the current identical harness and environment fixes. The optimized code applies
2x2 interaction blocks directly and scatters the UDT column permutation. The
rejected fused-rank-two candidate is absent from the final implementation.

The table below reports medians of **three sequential interleaved repetitions**.
`pairs` counts complete right+left sweeps. `speedup` is baseline/optimized median
sweep time; `paired range` shows the minimum and maximum ratio across matched
repetitions. `total speedup` includes initialization and one warmup sweep pair,
as well as the timed pairs. Validation and output serialization are excluded.
The comparison covers physical sites >100, not merely matrices with >100 rows.

See [the workflow and limitations](README.md), [raw final CSV](results/final.csv),
[runner log](results/benchmark.log), and [regression log](results/tests.log).

| Model | sites | matrix dim | Ltau | pairs | threads | seed | baseline sweep s | optimized sweep s | speedup | paired range | total speedup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| majorana | 128 | 128 | 12 | 8 | 1 | 42 | 1.252 | 0.828 | 1.51x | 1.49-1.56x | 1.41x |
| majorana | 288 | 288 | 12 | 8 | 1 | 42 | 13.757 | 8.581 | 1.60x | 1.60-1.62x | 1.50x |
| honeycomb | 128 | 256 | 12 | 8 | 1 | 42 | 9.682 | 5.950 | 1.63x | 1.60-1.65x | 1.52x |
| majorana | 512 | 512 | 12 | 4 | 1 | 42 | 41.683 | 27.460 | 1.52x | 1.52-1.53x | 1.41x |
| square | 144 | 288 | 12 | 4 | 1 | 42 | 8.861 | 5.770 | 1.54x | 1.50-1.59x | 1.42x |
| chain0 | 128 | 256 | 12 | 4 | 1 | 42 | 3.587 | 2.309 | 1.55x | 1.53-1.57x | 1.38x |
| chain1 | 128 | 256 | 12 | 4 | 1 | 42 | 3.603 | 2.365 | 1.52x | 1.50-1.58x | 1.39x |
| majorana | 288 | 288 | 12 | 4 | 4 | 923 | 2.792 | 2.039 | 1.37x | 1.32-1.44x | 1.27x |
| majorana | 128 | 128 | 40 | 4 | 1 | 123 | 2.062 | 1.422 | 1.45x | 1.43-1.45x | 1.36x |

Maximum raw signed-Pfaffian discrepancy: 6.191323814910013e-13
Maximum rebuilt Green relative error: 0.0


All 54 benchmark runs completed successfully. All 27 paired comparisons matched
auxiliary-field and Green-matrix checksums and complex signed phases within the
runner's tolerances. Median sweep speedups across the nine cases range from
1.37x to 1.63x; the corresponding measured total-runtime
speedups range from 1.27x to 1.52x.

The original baseline passed all 14 tests. The retained version passed all 18
regressions, including added dense-oracle, signed-update, and UDT tests, plus the
existing independent four-mode Fock/Pfaffian validator. The rejected fused
candidate passed 19 tests but lost on performance. Test logs are included under
[results](results/tests.log). Original executable smoke runs for Majorana
honeycomb and both chain HS channels are recorded in [smoke.log](results/smoke.log).

[Compiler flags, source/binary hashes and reference checks](results/provenance.txt)
record which code was measured. These are numerical and timing observations at
the stated parameters; they do not prove sign freedom, establish accuracy at
arbitrarily large inverse temperature, or measure independent samples per second.
