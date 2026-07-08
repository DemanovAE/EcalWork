#!/bin/bash

#SBATCH -p nica
#SBATCH -J ECAL
#SBATCH -a 0-9                 # 10 chunks: TASK_ID = 0..9
#SBATCH -N 1
#SBATCH -o /scratch3/dflusova/ECalWork/log/qa_%A_%a.out
#SBATCH -e /scratch3/dflusova/ECalWork/log/qa_%A_%a.err
#SBATCH -x ncx111,ncx112,ncx113,ncx115,ncx117,ncx121,ncx127,ncx146,ncx169,ncx171,ncx172,ncx181,ncx185,ncx203,ncx207,ncx208,ncx211,ncx213,ncx215,ncx216,ncx217,ncx222,ncx223,ncx224,ncx225,ncx226,ncx227

# --- Environment -------------------------------------------------------------

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v24.09.24-1
source /lhep/users/dflusova/mpdroot_polarization/mpdroot/install/config/env.sh


export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# --- Paths -------------------------------------------------------------------

# Build directory where EcalWorkExec lives
export BUILD_DIR=/lhep/users/dflusova/EcalWork/EcalWork/build
export EXEC=${BUILD_DIR}/EcalWorkExec

# Input ROOT file (single TTree, same for all tasks)
export INPUT_DIR=/nica/mpd1/demanov/ecal_mpd
export INPUT_FILE=${INPUT_DIR}/run_rc1-hs4_133_basket5_1616162.root

# Scratch for outputs/logs
export OUT_BASE=/scratch3/dflusova/ECalWork
export OUT_DIR=${OUT_BASE}/root
export LOG_DIR=${OUT_BASE}/log

export OUT_FILE=${OUT_DIR}/run_rc1-hs4_133_basket5_1616162_part${TASK_ID}.root
export LOG=${LOG_DIR}/log_qa_run8_${JOB_ID}_${TASK_ID}.log

mkdir -p "${OUT_DIR}"
mkdir -p "${LOG_DIR}"

# --- Define total entries and chunking ---------------------------------------
# IMPORTANT: set TOTAL_ENTRIES to the true GetEntries() of TTree "events"
TOTAL_ENTRIES=584947682      # from your log

N_CHUNKS=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))
CHUNK_SIZE=$(( TOTAL_ENTRIES / N_CHUNKS ))

FIRST_ENTRY=$(( TASK_ID * CHUNK_SIZE ))
if [ "${TASK_ID}" -eq "$((N_CHUNKS-1))" ]; then
  LAST_ENTRY=${TOTAL_ENTRIES}
else
  LAST_ENTRY=$(( (TASK_ID + 1) * CHUNK_SIZE ))
fi

# --- Logging -----------------------------------------------------------------

echo "Node name: ${SLURMD_NODENAME}"                                   &>> "${LOG}"
echo "INFILE:    ${INPUT_FILE}"                                        &>> "${LOG}"
echo "Job Id:    ${JOB_ID}"                                            &>> "${LOG}"
echo "Task Id:   ${TASK_ID}"                                           &>> "${LOG}"
echo "OUTFILE:   ${OUT_FILE}"                                          &>> "${LOG}"
echo "Entries:   total=${TOTAL_ENTRIES}"                               &>> "${LOG}"
echo "This task: firstEntry=${FIRST_ENTRY} lastEntry=${LAST_ENTRY}"    &>> "${LOG}"
echo "EXEC:      ${EXEC}"                                              &>> "${LOG}"

# --- Run executable ----------------------------------------------------------

cd "${BUILD_DIR}"

"${EXEC}" "${INPUT_FILE}" "${OUT_FILE}" ${FIRST_ENTRY} ${LAST_ENTRY} &>> "${LOG}"

echo "Job is finished" &>> "${LOG}"