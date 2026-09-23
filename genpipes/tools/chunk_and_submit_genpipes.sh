#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

CHUNK_SIZE=20
MAX_QUEUE=500
SLEEP_TIME=120
SCHEDULER=slurm
RETRY=10
SCHEDULER_USER=$USER

usage() {
  echo
  echo "usage: $0 <GENPIPES SCRIPT> <OUTPUT FOLDER> [options]"
  echo "  Chunk a Genpipes script and submit all chunks to the scheduler."
  echo
  echo "   <GENPIPES SCRIPT>       A Genpipes output script"
  echo "   <OUTPUT FOLDER>         Folder where to store chunks"
  echo "   -n <chunk size>         Maximum number of jobs per chunk, default=$CHUNK_SIZE"
  echo "   -q <MAX QUEUE>          Maximum number of jobs in scheduler queue, default=$MAX_QUEUE"
  echo "   -s <SLEEP TIME>         Seconds to sleep when queue is full, default=$SLEEP_TIME"
  echo "   -S <SCHEDULER>          Scheduler (slurm or pbs), default=$SCHEDULER"
  echo "   -l <N>                  Retry N times on chunk submit error, default=$RETRY"
  echo "   -u <USER>               Scheduler user to monitor queue for, default=$SCHEDULER_USER"
}

while getopts "hn:q:s:S:l:u:" opt; do
  case $opt in
    n) CHUNK_SIZE=${OPTARG} ;;
    q) MAX_QUEUE=${OPTARG} ;;
    s) SLEEP_TIME=${OPTARG} ;;
    S) SCHEDULER=${OPTARG} ;;
    l) RETRY=${OPTARG} ;;
    u) SCHEDULER_USER=${OPTARG} ;;
    h) usage; exit 0 ;;
    \?) usage; exit 1 ;;
  esac
done

shift $((OPTIND-1))

if [ $# -lt 2 ]; then
  usage
  exit 1
fi

genpipes_script=$1
out_dir=$2

if [ ! -f "$genpipes_script" ]; then
  echo "ERROR: genpipes script not found: $genpipes_script" >&2
  exit 1
fi

echo "=== Chunking genpipes script ==="
"${SCRIPT_DIR}/chunk_genpipes.sh" -n "$CHUNK_SIZE" "$genpipes_script" "$out_dir" || exit 1

echo "=== Submitting chunks ==="
"${SCRIPT_DIR}/submit_genpipes.sh" \
  -n "$MAX_QUEUE" \
  -s "$SLEEP_TIME" \
  -S "$SCHEDULER" \
  -l "$RETRY" \
  -u "$SCHEDULER_USER" \
  "$out_dir"
