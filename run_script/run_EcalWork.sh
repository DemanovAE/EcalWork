#!/bin/bash

#SBATCH -p nica
#SBATCH -J ECAL
#SBATCH -a 0-9                 # 10 chunks: TASK_ID = 0..9
#SBATCH -N 1
#SBATCH -o /scratch3/dflusova/ECalWork/log/qa_%A_%a.out
#SBATCH -e /scratch3/dflusova/ECalWork/log/qa_%A_%a.err
#SBATCH -x ncx113,ncx114,ncx116,ncx122,ncx123,ncx124,ncx125,ncx126,ncx132,ncx138,ncx142,ncx143,ncx144,ncx145,ncx146,ncx147,ncx148,ncx150,ncx151,ncx152,ncx153,ncx154,ncx155,ncx156,ncx157,ncx158,ncx159,ncx160,ncx161,ncx162,ncx163,ncx164,ncx165,ncx166,ncx167,ncx168,ncx170,ncx174,ncx175,ncx176,ncx177,ncx180,ncx182,ncx184,ncx186,ncx187,ncx188,ncx206,ncx212,ncx214,ncx228

# --- Environment -------------------------------------------------------------

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v24.09.24-1
source /lhep/users/dflusova/mpdroot_polarization/mpdroot/install/config/env.sh

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# --- Paths -------------------------------------------------------------------

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

# --- Define total entries and chunking --------------------------------------

TOTAL_ENTRIES=584947682      # from your log

N_CHUNKS=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))
CHUNK_SIZE=$(( TOTAL_ENTRIES / N_CHUNKS ))

FIRST_ENTRY=$(( TASK_ID * CHUNK_SIZE ))
if [ "${TASK_ID}" -eq "$((N_CHUNKS-1))" ]; then
  LAST_ENTRY=${TOTAL_ENTRIES}
else
  LAST_ENTRY=$(( (TASK_ID + 1) * CHUNK_SIZE ))
fi

# --- Selection criteria (steering) ------------------------------------------

# 0 = longitudinal analysis, 1 = transverse analysis
TRANS_FLAG=1

# Use QA histos for long / trans (1 = yes, 0 = no)
USE_LONG_QA=0
USE_TRANS_QA=1

LONG_MAX1=3500
LONG_MAX2=1500
MIN3X3=0.75
MAX3X3=0.99
MAXDIFF5X5=0.25

# Transverse cuts (steering)
TRANS_AMP1=300
TRANS_AMP2=500
TRANS_MINLEN=5
TRANS_CONTAM=0.50   # reserved, not used for now

# --- Logging -----------------------------------------------------------------

echo "Node name:      ${SLURMD_NODENAME}"                                 &>> "${LOG}"
echo "INFILE:         ${INPUT_FILE}"                                      &>> "${LOG}"
echo "Job Id:         ${JOB_ID}"                                          &>> "${LOG}"
echo "Task Id:        ${TASK_ID}"                                         &>> "${LOG}"
echo "OUTFILE:        ${OUT_FILE}"                                        &>> "${LOG}"
echo "Entries total:  ${TOTAL_ENTRIES}"                                   &>> "${LOG}"
echo "Entry range:    first=${FIRST_ENTRY} last=${LAST_ENTRY}"           &>> "${LOG}"
echo "EXEC:           ${EXEC}"                                            &>> "${LOG}"
echo "TRANS_FLAG:     ${TRANS_FLAG}"                                      &>> "${LOG}"
echo "USE_LONG_QA:    ${USE_LONG_QA}"                                     &>> "${LOG}"
echo "USE_TRANS_QA:   ${USE_TRANS_QA}"                                    &>> "${LOG}"
echo "LONG_MAX1:      ${LONG_MAX1}"                                       &>> "${LOG}"
echo "LONG_MAX2:      ${LONG_MAX2}"                                       &>> "${LOG}"
echo "MIN3X3:         ${MIN3X3}"                                          &>> "${LOG}"
echo "MAX3X3:         ${MAX3X3}"                                          &>> "${LOG}"
echo "MAXDIFF5X5:     ${MAXDIFF5X5}"                                      &>> "${LOG}"

# --- Write human-readable analysis summary ---------------------------------

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
  if [ "${TRANS_FLAG}" -eq 1 ]; then
  echo "Selection cuts (transverse strips):"
  echo "  Channel amplitude threshold_1 (ADC):      > ${TRANS_AMP1}"
  echo "  Strip neighbor threshold_2 (ADC):         > ${TRANS_AMP2}"
  echo "  Minimal strip length (cells):             >= ${TRANS_MINLEN}"
  echo "  Contamination fraction (reserved):        ${TRANS_CONTAM}"
  echo ""
  fi
  echo "Notes:"
  echo "  - \"hottest cell\" = cell with maximum integral in event."
  echo "  - \"second hottest\" = largest neighbour in 3x3 window."
  echo "  - Longitudinal QA histos live in directory 'longitudinal_QA' in the ROOT file."
  echo "  - Transverse QA histos live in directory 'transverse_QA' in the ROOT file."
} > "${SUMMARY_FILE}"

# --- Run executable ---------------------------------------------------------

cd "${BUILD_DIR}"

"${EXEC}" "${INPUT_FILE}" "${OUT_FILE}" ${FIRST_ENTRY} ${LAST_ENTRY} \
          ${TRANS_FLAG} ${USE_LONG_QA} ${USE_TRANS_QA} \
          ${LONG_MAX1} ${LONG_MAX2} ${MIN3X3} ${MAX3X3} ${MAXDIFF5X5} \
          ${TRANS_AMP1} ${TRANS_AMP2} ${TRANS_MINLEN} ${TRANS_CONTAM} &>> "${LOG}"

echo "Job is finished" &>> "${LOG}"