
# Functions

SLEEP_TIME=120
MAX_QUEUE=500
SHEDULER_USER=$USER
SCHEDULER=slurm

usage (){

echo
echo "usage: $0 <CHUNK FOLDER>
  Control the number of jobs in scheduler queue, and resubmit jobs if then fail at submit time."
echo
echo "   <CHUNK FOLDER>          The output folder from the chunk_genpipes.sh script"
echo "   -n <MAX QUEUE>          Maximum number of job in slurm queue"
echo "                             default=$MAX_QUEUE"
echo "   -s <SLEEP TIME>         Number of second to sleep when queue is full default=$SLEEP_TIME"
echo "   -S <SCHEDULER>          Scheduler running on the cluster (slurm or pbs) default=$SCHEDULER"

}



get_n_jobs () {
  if [[ ${SCHEDULER} ==  'slurm' ]]; then
    echo $(squeue -u $SHEDULER_USER -h -t pending,running | wc -l)
  elif [[ ${SCHEDULER} ==  'pbs' ]]; then
    echo $(showq  -u $SHEDULER_USER  | grep $SHEDULER_USER  | wc -l )
  fi
}

cancel_jobs () {
  echo ""
  job_script=$1
  echo cancel all jobs from ${job_script%.sh}.out

  if [[ ${SCHEDULER} ==  'slurm' ]]; then
    scancel $(cat ${job_script%.sh}.out | awk -F'=' '{print $2}')
  elif [[ ${SCHEDULER} ==  'pbs' ]]; then
    qdel $(cat ${job_script%.sh}.out | awk -F'=' '{print $2}')
  fi
  rm ${job_script%.sh}.out
}

cancel_trap () {
    cancel_jobs "$@"
    exit 0
}

submit () {
  job_script=${1}
  while true; do
    # clean cancel if there is an interruption
    trap "cancel_trap ${job_script%.sh}.out" EXIT
    bash ${job_script}
    ret_code=$?
    if [ ${ret_code} -eq 0 ]; then
      trap - SIGTERM
      touch ${job_script%.sh}.done
      echo ${job_script} was sucssfully submitted
      break
    else
      echo cancel all jobs from ${job_script%.sh}.out
      cancel_jobs ${job_script%.sh}.out
    fi
  done
}


#  Script

while getopts "hn:u:s:S:" opt; do
  case $opt in
    u)
      SHEDULER_USER=${OPTARG}
    ;;
    s)
      SLEEP_TIME=${OPTARG}
    ;;
    S)
      SCHEDULER=${OPTARG}
      if [[ ${SCHEDULER} != 'slurm'  && ${SCHEDULER} != 'pbs' ]] ;then
        echo only slurm and pbs scheduler are supported
        usage
        exit 1
      fi

    ;;
    n)
      MAX_QUEUE=${OPTARG}
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
chunk_folder=$1

if [ ! -d  ${chunk_folder} ]; then
  echo ${chunk_folder} does not exist
  exit 1
fi
# sourcing to get the value of CHUNK_SIZE
source ${chunk_folder}/header.sh





mkdir $chunk_folder/.lockdir 2>/dev/null
ret_code=$?
if [[ $ret_code -ne 0 ]] ; then
  echo it seems that another $0 process is runnning, stop it
  echo before restating
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
      if [[ $((MAX_QUEUE-curent_n_jobs)) -gt $CHUNK_SIZE ]]; then
       submit ${sh_script}
       break
      else
        echo to many jobs, sleeping for $SLEEP_TIME sec
        sleep $SLEEP_TIME
      fi
    done

  fi

done
