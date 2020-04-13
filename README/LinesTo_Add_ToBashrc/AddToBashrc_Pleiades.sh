# Source global definitions
if [ -f /etc/bash_profile ]; then
        . /etc/bash_profile
fi

export WORKSPACE_PATH='/u/sventuri/WORKSPACE/'
export LIB_PATH='/u/sventuri/libraries/'
export PATH=$HOME:$PATH

export MODULEPATH=/usr/share/modules/modulefiles:/nasa/modulefiles/sles12/:/u/sventuri/libraries/modules-3.2.10/usr/local/Modules/3.2.10/moduleapps:/nasa/modulefiles-sles11s
type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh


#LAPACK
export LAPACK_LIBS=${LIB_PATH}'/LAPACK/install/lib64/liblapack.so.3'
export LD_LIBRARY_PATH=$LAPACK_LIBS:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LAPACK_LIBS:$DYLD_LIBRARY_PATH

#BLAS
export BLAS_LIBS=${LIB_PATH}'/LAPACK/install/lib64/libblas.so.3'
export LD_LIBRARY_PATH=$BLAS_LIBS:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$BLAS_LIBS:$DYLD_LIBRARY_PATH

# #OpenBLAS
# export OBLAS_INC='/u/sventuri/libraries/OpenBLAS/install/OpenBLAS-1.13-intel-16.0.2/include'
# export OBLAS_LIB='/u/sventuri/libraries/OpenBLAS/install/OpenBLAS-1.13-intel-16.0.2/lib64'
# export LD_LIBRARY_PATH=$OBLAS_LIB:$LD_LIBRARY_PATH


########################################################## COMPILERS #############################################################
#INTEL
function LOAD_INTEL {
export FORTRAN_COMPILER=ifort
export COMPILER_VER=$INTEL_VER
export FC=ifort
export F77=ifort
export CMAKE_Fortran_COMPILER=ifort
export CC=icc
export CXX=icc
}

function INTEL_NEW_LOAD {
#module purge
#module load texlive/2016
#module load comp-intel/2016.2.181
module load comp-intel/2018.3.222
#module load python3/Intel_Python_3.6_2018.3.222
export OMP_NUM_THREADS=10
#export INTEL_VER='intel-16.0.2'
export INTEL_VER='intel-18.0.3'
LOAD_INTEL
}
##################################################################################################################################


################################################################ CoarseAIR ##############################################################
function COARSEAIR_UPDATE {
INTEL_NEW_LOAD
cd $WORKSPACE_PATH/CoarseAIR/coarseair/
source ./scripts/building/BuildingInstalling.sh
}
###################################################################################################################################
