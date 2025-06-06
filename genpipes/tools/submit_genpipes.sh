#!/bin/bash
# Functions

SLEEP_TIME=120
MAX_QUEUE=500
SCHEDULER_USER=$USER
SCHEDULER=slurm
export SUBMIT_RETCODE=1
export RETRY=10
type squeue > /dev/null 2>&1 || SCHEDULER=pbs

usage (){

echo
echo "usage: $0 <CHUNK FOLDER>
  Control the number of jobs in scheduler queue, and resubmit jobs if then fail at submit time."
echo
echo "   <CHUNK FOLDER>          The output folder from the chunk_genpipes.sh script"
echo "   -n <MAX QUEUE>          Maximum number of job in slurm queue, default=$MAX_QUEUE"
echo "   -s <SLEEP TIME>         Number of second to sleep when queue is full default=$SLEEP_TIME"
echo "   -S <SCHEDULER>          Scheduler running on the cluster (slurm or pbs) default=$SCHEDULER"
echo "   -l <N>                  Will retry N time to resubmit a chunk if error occurs default=$RETRY"

}

get_n_jobs () {
  if [[ ${SCHEDULER} ==  'slurm' ]]; then
    echo "$(squeue -u $SCHEDULER_USER -h -t pending,running | wc -l)"
  elif [[ ${SCHEDULER} ==  'pbs' ]]; then
    echo "$(showq  -u $SCHEDULER_USER  | grep $SCHEDULER_USER  | wc -l )"
  fi
}

cancel_jobs () {
  echo ""
  job_list=$1
  echo "$job_list"
  echo "cancel all jobs from ${job_list}"
  if [[ ${SCHEDULER} ==  'slurm' ]]; then
    scancel $(cat ${job_list} | awk -F'=' '{print $2}')
  elif [[ ${SCHEDULER} ==  'pbs' ]]; then
    qdel $(cat ${job_list} | awk -F'=' '{print $2}')
  fi
  rm ${job_list}  2>/dev/null
  echo "canceled all jobs from ${job_list}"
}

cancel_trap () {
    cancel_jobs "$@" 2>/dev/null
    rm -rf $chunk_folder/.lockdir
    exit 1
}

submit () {
  echo "submitting $1"
  job_script=${1}
  job_list=${job_script%.sh}.out
  for ((N=1;N<=RETRY;N++)); do
    # clean cancel if there is an interruption
    trap "echo cleanup; cancel_trap ${job_list}" EXIT
    bash ${job_script} 2> ${job_script%.sh}.err
    ret_code=$?
    if [ ${ret_code} -eq 0 ]; then
      trap - SIGTERM
      touch ${job_list}
      echo "${job_script} was sucessfully submitted"
      SUBMIT_RETCODE=0
      break
    else
      SUBMIT_RETCODE=1
      echo "error in submits"
      cancel_jobs ${job_list}
      sleep 1
      echo "resubmitting"
    fi
  done
  if [[ ${SUBMIT_RETCODE} -eq 1 ]]; then
    echo "could not complete submit after $RETRY retry"
    echo "Error log in ${job_script%.sh}.err:"
    cat "${job_script%.sh}.err"
    return 1
  fi
}


#  Script

while getopts "hl:n:u:s:S:" opt; do
  case $opt in
    u)
      SCHEDULER_USER=${OPTARG}
    ;;
    s)
      SLEEP_TIME=${OPTARG}
    ;;
    S)
      SCHEDULER=${OPTARG}
      if [[ ${SCHEDULER} != 'slurm'  && ${SCHEDULER} != 'pbs' ]] ;then
        echo "only slurm and pbs scheduler are supported"
        usage
        exit 1
      fi

    ;;
    n)
      MAX_QUEUE=${OPTARG}
    ;;
    l)
      RETRY=${OPTARG}
    ;;
    h)
      usage
      exit 0
    ;;
    \?)
      usage
      exit 1
    ;;
  esac
done

shift $((OPTIND-1))

if [ $# -lt 1 ]; then
  usage
  exit 1
fi
chunk_folder=$(realpath "$1")

if [ ! -d  ${chunk_folder} ]; then
  echo "${chunk_folder} does not exist"
  exit 1
fi
# sourcing to get the value of GENPIPES_CHUNK_SIZE
source ${chunk_folder}/header.sh
set +e

mkdir ${chunk_folder}/.lockdir 2>/dev/null
ret_code=$?
if [[ $ret_code -ne 0 ]] ; then
  echo "it seems that another $0 process is runnning"
  echo "If you are sure that no other process in running, run 'rm -r ${chunk_folder}/.lockdir'"
  echo "and restart $0"
  exit 1
else
  trap "rm -rf $chunk_folder/.lockdir" EXIT
fi

all_sh=($(ls ${chunk_folder}/chunk*sh|sort -V))
all_done=($chunk_folder/chunk*done)

for sh_script in "${all_sh[@]}"; do
  done_script=${sh_script%.sh}.done
  if [ ! -f $done_script ]; then
    while true ; do
      curent_n_jobs=$(get_n_jobs)
      if [[ $((MAX_QUEUE-curent_n_jobs)) -gt $GENPIPES_CHUNK_SIZE ]]; then
        submit ${sh_script}
        ret_code=$?
        trap "rm -rf $chunk_folder/.lockdir" EXIT
        if [[ $ret_code -ne 0 ]]; then
          # one of the chunks has failed, stop the process here
          exit $ret_code
        fi
        touch ${done_script}
        break
      else
          echo "too many jobs, sleeping for $SLEEP_TIME sec"
          sleep $SLEEP_TIME
      fi
    done
  fi
done

if [[ ${SUBMIT_RETCODE} -eq 0 ]]; then
  echo "All done, uploading usage statistics"
  bash ${chunk_folder}/wget_call.sh
else
  exit 1
fi
