#!/usr/bin/env bash

# FUNCTIONS

usage (){

echo
echo "usage: $0 <GENPIPES SCRIPT> <OUTPUT FOLDER>
  Chunk Genpipes submitting script so there is njobs in them"
echo
echo "   <GENPIPES SCRIPT>       A Genpipes output script"
echo "   <OUTPUT FOLDER>         Folder where to store chunks"
echo "   -n <chunk size>         Maximum number of job in chunk, default=20"

}

load_previous_submit_id (){
  current_file=$1
cat << 'EOF' >> "${current_file}"
for file in $(ls ${SCRIPTPATH}/*out | sort -n ); do
source $file
done
EOF
}

start_new_chunk () {
  header=$(basename "$2")
  cat << EOF > "$1"
#!/usr/bin/env bash
SCRIPTPATH="\$( cd "\$(dirname "\$0")" >/dev/null 2>&1 ; pwd -P )"
echo \${SCRIPTPATH}/${header}
source \${SCRIPTPATH}/${header}
STEP=$3
mkdir -p \$JOB_OUTPUT_DIR/\$STEP
EOF
}

# SCRIPT

max_chunk=20
while getopts "hn:" opt; do
  case $opt in
    n)
      max_chunk=${OPTARG}
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

if [ $# -lt 2 ]; then
  usage
  exit 1
fi
genpipes_in=$1
out_dir=$2

rm -rf "${out_dir}"/chunk_* 2>/dev/null
mkdir -p "${out_dir}"
header=${out_dir}/header.sh
echo '# header for all chunks' > "${header}"
while read -r line ; do
    if [[ $line =~ ^STEP=.*$ ]]; then
      STEP=${line#STEP=}
      break
    elif [[ $line =~ ^[^[:space:]]*=[^[:space:]]  ]]; then
      line="export ${line}"
    fi

    echo "$line" >> "${header}"
done < "${genpipes_in}"
echo "export GENPIPES_CHUNK_SIZE=${max_chunk}" >> "${header}"

chunk=1
out_file=/dev/null
# create chunks
while IFS= read -r line ; do
    echo "$line" >> "$out_file"

    # Literal comparison against the string as it appears in the genpipes script
    # shellcheck disable=SC2016
    if [ "$line" == 'cd $OUTPUT_DIR' ]; then
      out_file=${out_dir}/chunk_${chunk}.sh
      start_new_chunk "${out_file}" "${header}" "${STEP}"

    elif [[ $line =~ ^STEP=.*$ ]]; then
      STEP=${line#STEP=}

    elif [[ $line =~ \#.JOB:.* ]]; then
      nb_job=$((nb_job+1))
      if [[ ${nb_job} -ge ${max_chunk} ]]; then
          nb_job=0  # reset counter
          echo '# END' >> "$out_file"
          chunk=$((chunk+1))
          out_file=${out_dir}/chunk_${chunk}.sh
          start_new_chunk "${out_file}" "${header}" "${STEP}"
          load_previous_submit_id "${out_file}"

          echo "$line" >> "${out_file}"
      fi

    elif [[ $line =~ echo.*\ \>\>..JOB_LIST ]]; then
      job_list_chunk=chunk_${chunk}.out
      # shellcheck disable=SC2001  # backreference \1 requires sed; ${var//} cannot replicate this
      job_id_export_line=$(echo "$line" | sed 's/echo\s"$\(.*\)\s$JOB_NAME.*/echo "export \1=$\1 "/g')
      echo "$job_id_export_line >> \${SCRIPTPATH}/${job_list_chunk}" >> "${out_file}"
    fi
done < "${genpipes_in}"

# Exclude the wget call from the last chunk; submit_genpipes.sh will run it after all chunks are submitted
echo "removing wget call from $out_file"
split_on="Call home with pipeline statistics"
if ! csplit -sf "${out_dir}/wget_call" -n 1 "$out_file" /"${split_on}"/-1; then
  echo "ERROR: telemetry pattern '${split_on}' not found in ${out_file}" >&2
  exit 1
fi
mv "${out_dir}/wget_call0" "${out_file}"
mv "${out_dir}/wget_call1" "${out_dir}/wget_call.sh"
