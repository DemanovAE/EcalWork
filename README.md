# MPD Data Converter / ECalWork

Библиотека и утилиты для:
1. конвертации бинарных файлов ECAL MPD (`.data`) в ROOT-файлы с TTree `events`;
2. анализа этих ROOT-файлов (продольный / поперечный анализ ECAL) локально и на кластере (SLURM job array), с управляемыми критериями отбора через командную строку / bash-скрипты.

Library and tools for:
1. converting MPD ECAL binary `.data` files to ROOT files with an `events` TTree;
2. running longitudinal / transverse ECAL analysis on these ROOT files, locally or on a SLURM cluster (job arrays), with selection cuts controlled from the shell.

---

## 1. Build

```bash
git clone https://github.com/DemanovAE/EcalWork.git
cd EcalWork

mkdir build && cd build
cmake ..
make -j4 EcalWorkExec
```

After building you get:
- core data conversion library: `build/libMPDDataConverter.so` (or `.dylib` on macOS);
- analysis executable: `build/EcalWorkExec`.

---

## 2. Step 1: `.data` → ROOT (`ConvertToRoot`)

Conversion from raw `.data` to a ROOT file with `events` TTree is done by:

```cpp
void ConvertToRoot(int targetEvent = -1,
                   std::vector<int> targetEvents = {},
                   std::string inputData  = "../run_rc-hs1_088.data",
                   std::string outputData = "out_all.root");
```

Typical interactive usage in ROOT:

```bash
root -l EcalWork.cpp
```

Inside ROOT:

```cpp
.L EcalWork.cpp+

ConvertToRoot(-1, {},
              "/path/to/run_rc-hs1_088.data",
              "run_rc-hs1_088.root");
```

What `ConvertToRoot` does:
- reads the raw `.data` file via `MpdDataConverter`;
- applies pedestal and waveform processing (see `EcalWork.cpp`);
- produces a ROOT file `outputData` with an `events` TTree (branches `eventNum`, `integral`, `amplitude`, `channelNums`, `channel_Z`, `channel_Phi`);
- can optionally save waveform histograms / pictures under `pict/`.

---

## 3. Step 2: ROOT TTree analysis (`EcalWork` / `EcalWorkExec`)

The analysis logic lives in `EcalWork(...)` and is wrapped by `EcalWorkExec`.

### 3.1. Analysis configuration (EcalConfig)

Selection criteria are passed via a configuration struct:

```cpp
struct EcalConfig {
  bool  transAnalysis;              // 0 = longitudinal, 1 = transverse
  bool  useLongQA;                  // fill longitudinal QA histos
  bool  useTransQA;                 // fill transverse QA histos
  int   long_max_1st_Integral;      // cut on 1st hottest cell integral
  int   long_max_2nd_Integral;      // cut on 2nd hottest cell (3x3 neighbour)
  float long_min_3x3_ratio;         // min allowed max/(max+2nd) in 3x3
  float long_max_3x3_ratio;         // max allowed max/(max+2nd)
  float long_max_diffusivity_5x5;   // max (5x5−3x3)/5x5
};
```

This config is passed into:

```cpp
void EcalWork(std::string inputDataTree,
              std::string outputData,
              Long64_t firstEntry,
              Long64_t lastEntry,
              const EcalConfig &cfg);
```

Inside `EcalWork`:
- `cfg.transAnalysis` chooses between `LongAnalysis` and `TransverseAnalysis`;
- `cfg.useLongQA`, `cfg.useTransQA` control booking/writing of QA histograms:
  - `longitudinal_QA` directory for long selection;
  - `transverse_QA` directory for transverse selection;
- the long‑analysis cuts use the numeric fields of `cfg` (no hard‑coded globals).

### 3.2. Command-line interface: `EcalWorkExec`

`EcalWorkExec` is a thin wrapper that parses CLI arguments and builds `EcalConfig`:

```bash
./EcalWorkExec \
  <input.root> <output.root> \
  <firstEntry> <lastEntry> \
  <transFlag> <useLongQA> <useTransQA> \
  <longMax1> <longMax2> <min3x3> <max3x3> <maxDiff5x5>
```

Arguments:

- `<input.root>` – input file with an `events` TTree.
- `<output.root>` – output analysis file.
- `<firstEntry>` `<lastEntry>` – entry indices, processed range `[firstEntry, lastEntry)` (clamped to `GetEntries()`).
- `<transFlag>` – `0` = longitudinal selection, `1` = transverse selection.
- `<useLongQA>`, `<useTransQA>` – `1` to book/fill/write QA histograms, `0` to disable.
- `<longMax1>` – minimum integral of the hottest cell in the event.
- `<longMax2>` – minimum integral of the second hottest cell (within 3×3 around the max).
- `<min3x3>`, `<max3x3>` – allowed window for ratio `max / (max + 2nd)` in 3×3.
- `<maxDiff5x5>` – upper cut on diffusivity `(5x5−3x3)/5x5`.

Example: small longitudinal run with QA enabled and tuned cuts:

```bash
cd build

./EcalWorkExec \
  /nica/mpd1/demanov/ecal_mpd/run_rc1-hs4_133_basket5_1616162.root \
  /scratch3/dflusova/ECalWork/root/test_long.root \
  0 100000 \
  0 1 0 3500 1500 0.75 0.99 0.25
```

The output contains:
- summary histograms (event, channel, amplitude, integral);
- `target_channels` directory with per‑channel histos:
  - longitudinal mode: “core energy” (sum of hottest + strongest neighbour);
  - transverse mode: “target integral” of the selected strip cell;
- optional QA directories: `longitudinal_QA` and/or `transverse_QA`.

---

## 4. Physics interpretation of longitudinal cuts

In longitudinal mode the selection is designed to pick events where a single minimum-ionizing particle (muon) deposits energy in a well‑defined core cell, with limited spatial spread.

Conceptually:

1. **Hottest cell (1st max integral)**  
   - This is the cell with the largest signal integral in the event.  
   - Cut `long_max_1st_Integral` ensures the core signal is above noise / pedestal and consistent with at least ~1 MIP (exact MIP scale depends on calibration).

2. **Second hottest cell (3×3 neighbour)**  
   - Among the 3×3 window around the hottest cell, the largest neighbour is the “second hottest”.  
   - Cut `long_max_2nd_Integral` removes events where the neighbour is too small (e.g. very asymmetric energy sharing or noise).

3. **3×3 ratio: `max / (max + 2nd)`**  
   - Measures how much of the 3×3 core is carried by the central hottest cell.  
   - High ratio → very localized core (good MIP candidate).  
   - Low ratio → energy split between cells (possible shower or multi‑hit).  
   - Cuts `[long_min_3x3_ratio, long_max_3x3_ratio)` define an acceptable window.

4. **5×5 diffusivity: `(5x5 − 3x3) / 5x5`**  
   - Compares energy in a 5×5 window to energy in the inner 3×3.  
   - Small value → most energy is inside 3×3 → narrow shower, likely a single track.  
   - Large value → a lot of energy in outer ring → diffuse or multi‑particle event.  
   - Cut `long_max_diffusivity_5x5` rejects events with too much peripheral energy.

The longitudinal “core energy” stored per channel is the sum of the hottest cell and its strongest neighbour in 3×3. This is the main observable for muon studies and calibration scans.

## Transverse selection steering

Transverse selection now uses the same external steering pattern as
longitudinal selection.

**Config fields (struct `EcalConfig`):**

- `trans_amp_thr1`  – basic amplitude threshold_1 (ADC) for entering channels
- `trans_amp_thr2`  – stronger threshold_2 (ADC) used for strip finding
- `trans_min_strip_len` – minimal number of consecutive cells in a strip
- `trans_contam_frac`   – reserved fraction (0–1) for future contamination logic

These parameters are:

1. Initialised with defaults in `main()`.
2. Overridable from command-line arguments.
3. In practice, set from the SLURM steering script.

Inside `TransverseAnalysis`:

- All previous hard-coded numbers (100, 500, 5, 0.20) are replaced by
  the corresponding `cfg.trans_*` fields.
- The physics logic is unchanged; only the source of the thresholds moved
  from literals to the config struct.
---

## 5. Cluster run (SLURM job array)

For large ROOT files, use the SLURM script in `run_script`, e.g. `run_EcalWork.sh`.

### 5.1. SLURM script structure (example)

Key parts of the script:

```bash
#!/bin/bash

#SBATCH -p nica
#SBATCH -J ECAL
#SBATCH -a 0-9                 # 10 chunks: TASK_ID = 0..9
#SBATCH -N 1
#SBATCH -o /scratch3/dflusova/ECalWork/log/qa_%A_%a.out
#SBATCH -e /scratch3/dflusova/ECalWork/log/qa_%A_%a.err
#SBATCH -x ncx111,ncx112,...

# Environment
source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v24.09.24-1
source /lhep/users/dflusova/mpdroot_polarization/mpdroot/install/config/env.sh

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# Paths
export BUILD_DIR=/lhep/users/dflusova/EcalWork/EcalWork/build
export EXEC=${BUILD_DIR}/EcalWorkExec

export INPUT_DIR=/nica/mpd1/demanov/ecal_mpd
export INPUT_FILE=${INPUT_DIR}/run_rc1-hs4_133_basket5_1616162.root

export OUT_BASE=/scratch3/dflusova/ECalWork
export OUT_DIR=${OUT_BASE}/root
export LOG_DIR=${OUT_BASE}/log

export OUT_FILE=${OUT_DIR}/run_rc1-hs4_133_basket5_1616162_part${TASK_ID}.root
export LOG=${LOG_DIR}/log_qa_run8_${JOB_ID}_${TASK_ID}.log

mkdir -p "${OUT_DIR}"
mkdir -p "${LOG_DIR}"
```

Chunking logic:

```bash
TOTAL_ENTRIES=584947682    # events->GetEntries() for this input file

N_CHUNKS=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))
CHUNK_SIZE=$(( TOTAL_ENTRIES / N_CHUNKS ))

FIRST_ENTRY=$(( TASK_ID * CHUNK_SIZE ))
if [ "${TASK_ID}" -eq "$((N_CHUNKS-1))" ]; then
  LAST_ENTRY=${TOTAL_ENTRIES}
else
  LAST_ENTRY=$(( (TASK_ID + 1) * CHUNK_SIZE ))
fi
```

Selection steering (example):

```bash
TRANS_FLAG=0       # 0 = longitudinal, 1 = transverse
USE_LONG_QA=1
USE_TRANS_QA=0

LONG_MAX1=3500
LONG_MAX2=1500
MIN3X3=0.75
MAX3X3=0.99
MAXDIFF5X5=0.25
```

Run:

```bash
cd "${BUILD_DIR}"

"${EXEC}" "${INPUT_FILE}" "${OUT_FILE}" ${FIRST_ENTRY} ${LAST_ENTRY} \
          ${TRANS_FLAG} ${USE_LONG_QA} ${USE_TRANS_QA} \
          ${LONG_MAX1} ${LONG_MAX2} ${MIN3X3} ${MAX3X3} ${MAXDIFF5X5} &>> "${LOG}"

echo "Job is finished" &>> "${LOG}"
```

Each array task (`TASK_ID = 0..9`) processes its own entry range and writes one ROOT and one log file.

### 5.2. Human-readable analysis summary per job

The script can also write a small text file describing the configuration, e.g.:

```bash
SUMMARY_FILE=${OUT_DIR}/technical_notes_about_jobs/analysis_summary_part${TASK_ID}.txt

{
  echo "=== ECal cosmic analysis summary ==="
  echo ""
  echo "Date/time:        $(date)"
  echo "Node:             ${SLURMD_NODENAME}"
  echo "Job ID:           ${JOB_ID}"
  echo "Task ID:          ${TASK_ID}"
  echo ""
  echo "Input file:       ${INPUT_FILE}"
  echo "Output file:      ${OUT_FILE}"
  echo ""
  echo "Total entries:    ${TOTAL_ENTRIES}"
  echo "Entry range:      [${FIRST_ENTRY}, ${LAST_ENTRY})"
  echo ""
  if [ "${TRANS_FLAG}" -eq 0 ]; then
    echo "Analysis mode:    Longitudinal selection"
  else
    echo "Analysis mode:    Transverse selection"
  fi
  echo ""
  echo "QA histograms:"
  echo "  Longitudinal QA: $( [ "${USE_LONG_QA}" -eq 1 ] && echo ENABLED || echo disabled )"
  echo "  Transverse QA:   $( [ "${USE_TRANS_QA}" -eq 1 ] && echo ENABLED || echo disabled )"
  echo ""
  echo "Selection cuts (longitudinal core if used):"
  echo "  1st hottest cell integral (ADC):           > ${LONG_MAX1}"
  echo "  2nd hottest cell integral in 3x3 (ADC):    > ${LONG_MAX2}"
  echo "  3x3 ratio = max / (max + 2nd):             [${MIN3X3}, ${MAX3X3})"
  echo "  Diffusivity = (5x5 - 3x3) / 5x5:           < ${MAXDIFF5X5}"
  echo ""
  echo "Notes:"
  echo "  - \"hottest cell\" = cell with maximum integral in event."
  echo "  - \"second hottest\" = largest neighbour in 3×3 window."
  echo "  - Longitudinal QA histos live in 'longitudinal_QA' in the ROOT file."
  echo "  - Transverse QA histos live in 'transverse_QA' in the ROOT file."
} > "${SUMMARY_FILE}"
```

---

## 6. Merging output ROOT files

After the array finishes, you have:

```text
run_rc1-hs4_133_basket5_1616162_part0.root
...
run_rc1-hs4_133_basket5_1616162_part9.root
```

You can merge histograms with `hadd`:

```bash
hadd -f run_rc1-hs4_133_basket5_1616162_all.root \
     run_rc1-hs4_133_basket5_1616162_part*.root
```

or use a custom ROOT script if you need more control over directories.

---

## 7. Requirements

- ROOT 6
- C++17
- MPDROOT environment (e.g.):
  - `source /cvmfs/nica.jinr.ru/sw/os/login.sh latest`
  - `module add mpddev/v24.09.24-1`
  - `source /lhep/users/dflusova/mpdroot_polarization/mpdroot/install/config/env.sh`
- SLURM for cluster runs (if using job arrays).
