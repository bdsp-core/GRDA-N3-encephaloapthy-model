# Restoring states 49–106

States 49–106 were restored on branch `states-49-106`, based on `8412e02`.

## Source and scope

The source is the historical uncommitted `config.py`, `cell_classes.py`, `network_class.py`, and `bazh_net.py` in `/Users/GrahamMcLeod/GRDA_sleep_model/fink_sleep`. The historical config SHA256 is recorded in `preset_parity.json` and the committed-test fixture `tests/fixtures/historical_state_parameters.json`.

- Restored every definition for states 49–106, including defined but previously unrun states. State 4 remains reserved. States 107 onward remain outside this handoff.
- Exposed the historical route-specific AMPA and random-event amplitude controls, TC Ih conductance and TC T-current conductance. State 90 sets TC T-current conductance to zero. States 91–97 disable specified pathways by setting conductance to zero; connections/NetCons remain present.
- Restored optional TC hyperpolarizing current and diagnostic vectors, bias-specific output filenames, and protection against overwriting biased runs.
- Preserved main's optimized `CorticalFieldRecorder`, combined MPI gather, vectorized output writer, configurable threads, and existing sleep-state sweep switch.
- All `.mod` files are unchanged from `8412e02`: THREADSAFE and Rinf/Rtau RANGE changes are retained. No older `.mod` files were copied.
- Historical current/diagnostic vectors needed explicit owning point-process arguments to support main's newer multi-threaded runs. This is the only additional implementation adaptation to the restored interventions.

## How to reproduce the TC bias

The bias is independent of state ID. Merely choosing state 99, 103, or 106 does **not** apply it. Historical biased runs use:

- `FINK_TC_ICLAMP_UA_CM2=-0.30` (µA/cm²)
- `FINK_TC_ICLAMP_RAMP_MS=5000` (linear ramp from zero)
- `FINK_RECORD_TC_VOLTAGE=1` (required with nonzero bias)

Each TC cell has 29,000 µm² membrane area: −0.30 µA/cm² gives −0.087 nA per cell. The filename suffix is `_tcIm0p3_r5000ms` before `_nhost=N`.

After activating the environment and compiling with `nrnivmodl mod`, a launch using the existing sweep runner is:

```bash
FINK_TC_ICLAMP_UA_CM2=-0.30 \
FINK_TC_ICLAMP_RAMP_MS=5000 \
FINK_RECORD_TC_VOLTAGE=1 \
python run_sweep.py --states 99 103 106 --seeds 1-3 --duration-ms 120000 --ranks 8 --threads 1 --jobs 1
```

Choose seeds and MPI rank count from the supplied historical manifest for any specific file. This example uses the eight MPI ranks recorded in the historical manifest, with one simulation at a time. `FINK_TC_ICLAMP_ALLOW_OVERWRITE=1` explicitly permits replacement of an existing biased run.

## Verification

- Every `init_*` parameter of all 106 supported state IDs from 0 through 106 exactly matches the historical source (state 4 excluded). Fixture is generated from historical source, independently of restored code.
- Every original scalar configuration value in all 48 supported states 0–48 exactly matches the fresh-main source.
- `nrnivmodl mod` compiles successfully with NEURON 9.0.1.
- `FINK_RUN_NEURON_TESTS=1 python -m unittest discover -s tests`: **35 tests passed**. These include actual-cell tests for zero TC T-current, cortex–thalamus lesions preserving the internal TC–RE loop, and −0.087 nA delivered to all 100 TC cells with two NEURON threads.
- A complete 100-ms run of state 99, one MPI rank/two threads, passed. It produced 4,000 LFP and 4,000 vcort samples and diagnostics for eight TC cells at all 101 expected 1-ms timestamps. A second run is correctly blocked from overwriting it.
- A complete 100-ms state-106 run with **two MPI ranks and two threads per rank** also passed the same output checks.
- `git diff --check` passes.

The brief smoke runs used a 10-ms ramp so that the clamp reaches its target within 100 ms. They verify integration and output completeness, not the long-run dynamical phenotype. No long historical simulation was rerun, and byte-for-byte agreement of new runs with the historical uncommitted outputs is **not established**. Matching the MPI rank count alone is insufficient to promise identical output across code or runtime versions.

Saved validation logs and smoke traces are under `outputs/eli_handoff/validation` (outside the publication clone).
