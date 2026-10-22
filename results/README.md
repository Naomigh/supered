# Benchmark output

Run benchmarks into new subdirectories here. Generated results and builds are
excluded from Git. The original development workspace retains the large raw
inputs, per-window records, executable snapshots, and historical reports;
those archives are not required to build or use the library.

For a fresh clone, build `benchmark_single_pair` and `cache_control`, then use
`scripts/run_single_pair.py NEW_OUTPUT_DIRECTORY` and
`scripts/analyze_single_pair.py NEW_OUTPUT_DIRECTORY`. This benchmark uses
machine-specific Linux perf/RDPMC events and pins CPU 24; consult
`bench/SINGLE_PAIR_PROTOCOL.txt` before running it on another machine.

Legacy profiling and report scripts are included for provenance. Some require
the original development workspace and its earlier archived experiments.
