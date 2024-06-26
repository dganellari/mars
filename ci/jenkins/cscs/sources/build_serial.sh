set -e
set -x
set -o pipefail

BUILDBASE=${1:-${SCRATCH}/build/serial}
INSTALLBASE=${2:-${PROJECT}/install/serial}

SCRIPTBASE=$(dirname $(realpath $0))

unset CMAKE_DIR
unset TRILINOS_DIR

mkdir -p "${INSTALLBASE}"
rm -f "${INSTALLBASE}/environment"
#echo "source /opt/cray/pe/cdt/19.10/restore_system_defaults.sh >/dev/null 2>&1" > "${INSTALLBASE}/environment"
echo "module purge >/dev/null 2>&1" >> "${INSTALLBASE}/environment"
MODULES="modules cce daint-gpu PrgEnv-gnu slurm cray-mpich cray-hdf5-parallel cray-netcdf-hdf5parallel cray-libsci cray-python " #VTK"
for m in ${MODULES} ; do
  echo "module load ${m}" >> "${INSTALLBASE}/environment"
done
echo "export CRAYPE_LINK_TYPE=dynamic" >> "${INSTALLBASE}/environment"
echo "export CC=cc" >> "${INSTALLBASE}/environment"
echo "export CXX=CC" >> "${INSTALLBASE}/environment"
echo "export FC=ftn" >> "${INSTALLBASE}/environment"
echo "export F90=ftn" >> "${INSTALLBASE}/environment"
echo "export F77=ftn" >> "${INSTALLBASE}/environment"
echo "CXXFLAGS=\"\$(printf '%s\n' \"\${CXXFLAGS} -std=c++14\" | awk -v RS='[[:space:]]+' '!a[\$0]++{printf \"%s%s\", \$0, RT}')\"" >> "${INSTALLBASE}/environment"
echo "export CXXFLAGS" >> "${INSTALLBASE}/environment"

ENVBASE=/apps/daint/UES/anfink/cpu
export CMAKE_DIR="${ENVBASE}/cmake"
echo "export CMAKE_DIR=\"${CMAKE_DIR}\"" >> "${INSTALLBASE}"/environment
echo "export PATH=\"\${CMAKE_DIR}/bin:\${PATH}\"" >> "${INSTALLBASE}"/environment

source "${INSTALLBASE}/environment"

mkdir -p ${BUILDBASE}/mars
"${SCRIPTBASE}/mars.sh" "${BUILDBASE}/mars" "${INSTALLBASE}/mars" |& tee ${BUILDBASE}/mars/logfile
export MARS_DIR="${INSTALLBASE}/mars"
cp ${BUILDBASE}/mars/logfile ${INSTALLBASE}/mars/build.log
rm -Rf "${BUILDBASE}/mars"
echo "export MARS_DIR=\"${MARS_DIR}\"" >> "${INSTALLBASE}/environment"

echo "export LD_LIBRARY_PATH=\"\${TRILINOS_DIR}/lib:\$LD_LIBRARY_PATH\"" >> "${INSTALLBASE}/environment"
