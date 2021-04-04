#!/bin/bash -l
#SBATCH --job-name=mars_build_jenkins
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --partition=cscsci
#SBATCH --constraint=gpu
#SBATCH --error=/scratch/snx3000/gandanie/build/jenkins-mars%j.err
#SBATCH --output=/scratch/snx3000/gandanie/build/jenkins-mars%j.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

set -e
set -x

id
ls -lh
pwd -P
env
hostname
setfacl --version
DATE=$(date +"%Y%m%d%H%M")
SCRATCH=/scratch/snx3000/gandanie

function set_rights {
  echo "Setting rights"
  setfacl -R -m m::rwx ${BUILDDIR}
  setfacl -R -m m::rwx ${INSTALLDIR}
}
trap set_rights ERR

# Delete old builds. We need to loop over the results of find, instead of passing it to a `-exec`, otherwise the find command fails
# This command fails, since it tries to descend into {}.install after it has been deleted
# find $SCRATCH/build -mindepth 2 -maxdepth 2 -type d -regex '.*20.*\(Debug\|Release\)' -exec rm -Rf {} {}.install \;
for dir in $(find $SCRATCH/build -mindepth 2 -maxdepth 2 -type d -regex '.*20.*\(Debug\|Release\)') ; do
  rm -Rf "${dir}" "${dir}.install"
done

#for cpu_gpu in mars_cpu mars_gpu ; do
for cpu_gpu in cpu gpu; do
  for debug_release in Debug Release ; do
    BUILDDIR=${SCRATCH}/build/${cpu_gpu}/${DATE}_${BUILD_NUMBER}_${debug_release}
    INSTALLDIR=${BUILDDIR}.install

    built_new_version="false"
    if [[ "${debug_release}" == "Debug" ]]; then
      if [[ $(( ${BUILD_NUMBER} % 7 )) == 3 ]]; then
        ~gandanie/buildscripts/build_${cpu_gpu}_debug.sh ${BUILDDIR} ${INSTALLDIR}
        built_new_version="true"
      fi
    else
      ~gandanie/buildscripts/build_${cpu_gpu}.sh ${BUILDDIR} ${INSTALLDIR}
      built_new_version="true"
    fi

    find $(dirname ${BUILDDIR}) -maxdepth 1 -mindepth 1 -user jenkssl -name "*${debug_release}" -exec rm -Rf {} +
    DAYS="+1"
    [[ ${debug_release} == "Debug" ]] && DAYS="+7"
    find $(dirname ${INSTALLDIR}/) -maxdepth 1 -mindepth 1 -user jenkssl -mtime ${DAYS} -name "*${debug_release}.install" -exec rm -Rf {} +

    if [[ ${built_new_version} == "true" ]]; then
      LINK="$(dirname ${INSTALLDIR})/${cpu_gpu}_jenkins_${debug_release}"
      rm -f "${LINK}"
      ln -s "${INSTALLDIR}" "${LINK}"
    fi
  done
done
