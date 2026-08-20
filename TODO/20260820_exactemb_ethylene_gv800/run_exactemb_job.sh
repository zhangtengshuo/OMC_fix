#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

[[ $# -eq 3 ]] || { echo "Usage: $0 PROJECT RUN_DIR MOLCAS_MEM_MIB" >&2; exit 2; }
PROJECT=$1
RUN_DIR=$2
MOLCAS_MEM_MIB=$3
[[ "$(basename "$RUN_DIR")" == "$PROJECT" ]] || { echo "project/run-dir mismatch" >&2; exit 2; }
[[ -f "$RUN_DIR/$PROJECT.in" ]] || { echo "missing input: $RUN_DIR/$PROJECT.in" >&2; exit 2; }
[[ ! -e "$RUN_DIR/$PROJECT.out" ]] || { echo "refusing to overwrite output: $RUN_DIR/$PROJECT.out" >&2; exit 3; }

source /opt/profile.d/noci-openmolcas.source
export MOLCAS=/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11
export MOLCAS_LIBMOLCAS_SO="$MOLCAS/lib/libmolcas.so"
export PATH="$MOLCAS:$MOLCAS/bin:/opt/mamba/envs/omc/bin:$PATH"
export LD_LIBRARY_PATH="$MOLCAS/lib:/scratch/alpha/build-exactden-20260820_093910/wignernj-sys/lib:/opt/mamba/envs/omc/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export MOLCAS_NPROCS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"
export MOLCAS_MEM="$MOLCAS_MEM_MIB"
export MOLCAS_KEEP_WORKDIR=YES
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMPI_MCA_plm='^slurm'
export OMPI_MCA_ras='^slurm'
unset DISPLAY

EXPERIMENT_ROOT=/scratch/alpha/exactemb-ethylene-20260820_103419
export MOLCAS_SCRATCH_ROOT="$EXPERIMENT_ROOT/work/openmolcas"
export MOLCAS_WORKDIR="$MOLCAS_SCRATCH_ROOT/${PROJECT}.${SLURM_JOB_ID:?SLURM_JOB_ID is not set}"

cleanup() {
  local cleanup_rc=$?
  trap - EXIT
  if [[ -e "$MOLCAS_WORKDIR" ]]; then
    find "$MOLCAS_WORKDIR" -depth -delete || cleanup_rc=91
  fi
  exit "$cleanup_rc"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

mkdir -p "$MOLCAS_WORKDIR"
cd "$RUN_DIR"
{
  echo "Job ID: $SLURM_JOB_ID"
  echo "Host: $(hostname)"
  echo "Project: $PROJECT"
  echo "OpenMolcas: $MOLCAS"
  echo "MPI processes: $MOLCAS_NPROCS"
  echo "Memory per process MiB: $MOLCAS_MEM"
  echo "Task scratch: $MOLCAS_WORKDIR"
  date --iso-8601=seconds
} > "$PROJECT.provenance"

set +e
/usr/bin/time -f 'WALL_SECONDS=%e\nMAX_RSS_KB=%M' "$MOLCAS/pymolcas" "$PROJECT.in" > "$PROJECT.out" 2>&1
run_rc=$?
set -e
printf 'RUN_RC=%s\n' "$run_rc" > "$PROJECT.status"
exit "$run_rc"
