
# Functions

usage (){

echo
echo "usage: $0 <CHUNK FOLDER>
  Chunk Genpipes submiting sript so there is njobs in them"
echo
echo "   <Genpipes script>       A Genpipes output script"
echo "   -n <MAX QUEUE>          Maximum number of job in slurm queue"
echo "                             default=500"
echo "   -s <SLEEP TIME>          number of second to sleep when queue is full default= 120"

}

cancel_jobs () {
  echo ""
}

submit () {
  job_script=${1}
  while true; do
    bash ${job_script}
    ret_code=$?
    if [ ${ret_code} -eq 0 ]; then
      touch ${job_script%.sh}.done
      echo ${job_script} was sucssfully submitted
      break
    else
      echo cancel all jobs from ${job_script%.sh}.out
      scancel $(cat ${job_script%.sh}.out | awk -F'=' '{print $2}')
      rm ${job_script%.sh}.out
    fi
  done
}


#  Script

SLEEP_TIME=120
max_queue=500
slurm_user=$USER
while getopts "hn:u:s:" opt; do
  case $opt in
    u)
      slurm_user=${OPTARG}
    ;;
    s)
      SLEEP_TIME=${OPTARG}
    ;;
    n)
      max_queue=${OPTARG}
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
# all_out=$chunk_folder/*out
all_done=($chunk_folder/chunk*done)

for sh_script in ${all_sh[@]}; do
  done_script=${sh_script%.sh}.done
  echo ${done_script}
  if [ ! -f $done_script ]; then

    while true ; do
      curent_job_n=$(squeue -u $slurm_user -h -t pending,running | wc -l)
      if [[ $((max_queue-$curent_job_n)) -gt $CHUNK_SIZE ]]; then
       submit ${sh_script}
       break
      else
        echo to many jobs, sleeping for $SLEEP_TIME sec
        sleep $SLEEP_TIME
      fi
    done

  fi

done
