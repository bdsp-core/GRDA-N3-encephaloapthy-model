# Comparing this NEURON port against the original C++

The NEURON port produces an N2 oscillation at about 6.9 Hz, below the 7-15 Hz
the paper states for the model and well below the 11-16 Hz of human sleep
spindles. This directory reproduces the check of whether that comes from the
port or from the original model.

**Result: it comes from the original model.** Measuring the same quantity in
both implementations, the dominant N2 frequency is 6.62 Hz in the C++ and
6.88 Hz in the NEURON port. See [../docs/FINDINGS.md](../docs/FINDINGS.md).

## Getting and building the original

```bash
git clone https://github.com/bazhlab-ucsd/sleep-stage-transition cpp
cd cpp
```

Two changes are needed to build on macOS. Neither touches the model.

**1. Reentrant RNG.** The code uses glibc's `drand48_r`/`srand48_r`, which do
not exist on Darwin. `drand48_r_shim.h` in this directory implements the same
48-bit LCG that glibc uses, so the number stream is identical rather than
merely equivalent. Include it at build time; the model's own sources are left
untouched.

**2. A Linux-only diagnostic.** `print_used_mem()` in `io.cpp` reads
`/proc/self/status` and then calls `fgets` on the result without checking it.
On any non-Linux platform the `fopen` returns NULL and the program segfaults
before the simulation starts. Add a guard after the `char line[128];`
declaration:

```c
if (!proc) { printf("Resident RAM: unavailable on this platform\n"); return; }
```

Then build with a compiler that has OpenMP (Apple clang does not; `brew install gcc`):

```bash
g++-16 -O2 generate_network.cpp -o generate_network
./generate_network network.cfg > connection_info2

g++-16 -g -O3 -fopenmp -include ../drand48_r_shim.h \
    main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp -o mindcurrent -lm
```

## Selecting a single vigilance state

`main.cpp` sets `int one_state=0`, which runs the full ultradian cycle with
stage boundaries at `tmax/5`. The values 1-4 select a single state
(1 = awake, 2 = N2, 3 = N3, 4 = REM), but nothing assigns them: the
`one_state.h` file in the repository contains `4` and is never included by any
source file. To select a state without editing logic, make the variable read
the environment:

```c
int one_state = getenv("ONE_STATE") ? atoi(getenv("ONE_STATE")) : 0;
```

Then a static N2 run, matching this port's `FINK_STATE=1`:

```bash
sed -e 's/^tmax 350000/tmax 130000/' -e 's/^num_mp_threads 30/num_mp_threads 8/' \
    params.txt > params_n2.txt
mkdir -p out_n2
ONE_STATE=2 ./mindcurrent params_n2.txt out_n2 connection_info2
```

## Comparing

`out_n2/time_cx` holds the membrane voltage of all 500 cortical cells sampled
at 1 ms, one row per sample with time in the first column. The mean across
cells is the paper's own LFP definition. Compare it against this port's
`vcort_state1_seed1_nhost=8.txt`, which is the summed cortical voltage and so
the same quantity up to a scale factor that does not affect frequency:

```bash
../.venv/bin/python compare_spindle_frequency.py \
    --cpp-time-cx cpp/out_n2/time_cx \
    --port-vcort ../vcort_state1_seed1_nhost=8.txt
```

## Notes and caveats

- `params.txt` carries an upstream warning that changing the OpenMP thread
  count changes the output. Results are therefore not bitwise reproducible
  across thread counts; the comparison above used 8 threads.
- The C++ is slower than the NEURON port for this model: roughly 26x
  real-time on 8 threads, against roughly 7.5x real-time for the port on 8 MPI
  ranks. Most of the difference is that the C++ writes every cell's voltage
  every millisecond.
- The C++ writes per-cell voltage, which the NEURON port does not. That makes
  it useful for questions the port cannot answer directly, such as spatial
  synchrony across the cortical array.
