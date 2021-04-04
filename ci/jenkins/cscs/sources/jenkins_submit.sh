#!/bin/bash

JOB_ID=""

set -x

function cancel_job {
  echo "inside cancel_job"
  if [[ -n ${JOB_ID} ]]; then
    echo ${JOB_ID}
    JOB_ID="$(echo ${JOB_ID} | sed -e 's|[^0-9]*||')"
    echo $JOB_ID
    squeue -j ${JOB_ID} >& /dev/null && scancel ${JOB_ID}
  fi
}
trap cancel_job 15

JOB_ID="$(sbatch --wait $@)"
