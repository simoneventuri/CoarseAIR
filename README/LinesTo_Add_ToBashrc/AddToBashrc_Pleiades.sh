# Source global definitions
if [ -f /etc/bash_profile ]; then
        . /etc/bash_profile
fi

export WORKSPACE_PATH='/u/sventuri/WORKSPACE/'
export LIB_PATH='/u/sventuri/libraries/'

export MODULEPATH=/usr/share/modules/modulefiles:/nasa/modulefiles/sles12/:/u/sventuri/libraries/modules-3.2.10/usr/local/Modules/3.2.10/moduleapps:/nasa/modulefiles-sles11s
type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh

#alias ls='ls --icolor=tty'
alias lx='ls -lt'
alias vi='gvim'
alias lr='ls -lr'
alias lsd='ls -d */.'
alias cp='cp -i'
alias rm='rm -i'
alias cmake='/nasa/pkgsrc/sles12/2018Q3/bin/cmake'
export CLICOLOR=1
export LSCOLORS=GxFxCxDxBxegedabagaced
export FORT_FMT_RECL=4000
#module load texlive/2016
export PATH=$HOME:$PATH


#LAPACK
#LAPACK_LIBS='/usr/lib64/liblapack.so.3.5.0'
LAPACK_LIBS=${LIB_PATH}'/LAPACK/install/lib64/liblapack.so.3'
export LAPACK_LIBS
export LD_LIBRARY_PATH=$LAPACK_LIBS:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$LAPACK_LIBS:$DYLD_LIBRARY_PATH

#BLAS
#BLAS_LIBS='/usr/lib64/libblas.so.3.5.0'
BLAS_LIBS=${LIB_PATH}'/LAPACK/install/lib64/libblas.so.3'
export BLAS_LIBS
export LD_LIBRARY_PATH=$BLAS_LIBS:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$BLAS_LIBS:$DYLD_LIBRARY_PATH

#OpenBLAS
#OBLAS_INC='/u/sventuri/libraries/OpenBLAS/install/OpenBLAS-1.13-intel-16.0.2/include'
#export OBLAS_INC
#OBLAS_LIB='/u/sventuri/libraries/OpenBLAS/install/OpenBLAS-1.13-intel-16.0.2/lib64'
#export OBLAS_LIBS
#export LD_LIBRARY_PATH=$OBLAS_LIB:$LD_LIBRARY_PATH


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
module load comp-intel/2016.2.181
#module load cmake/2.8.12.1
export INTEL_VER='intel-16.0.2'
LOAD_INTEL
}
##################################################################################################################################


###################################################### GSL & FGSL ################################################################
#FGSL
function FGSL_for_COARSEAIR {
FGSL_DIR=${LIB_PATH}"/FGSL/install/fgsl-0.9.3-gsl-1.14-ifort-16.0.2/"
export FGSL_DIR
FGSL_INC=${LIB_PATH}"/FGSL/install/fgsl-0.9.3-gsl-1.14-ifort-16.0.2/include/ifort/"
export FGSL_INC
FGSL_LIB=${LIB_PATH}"/FGSL/install/fgsl-0.9.3-gsl-1.14-ifort-16.0.2/lib/libfgsl_ifort.a"
export FGSL_LIB
export LD_LIBRARY_PATH=$FGSL_LIB:$LD_LIBRARY_PATH
FGSL_LIB_FLAG='-lfgsl'
export FGSL_LIB_FLAG
}

function FGSL_for_SMUQ {
FGSL_DIR=${LIB_PATH}"/FGSL/install/fgsl-1.2.0-gsl-2.3-ifort-16.0.2/"
export FGSL_DIR
FGSL_INC=${LIB_PATH}"/FGSL/install/fgsl-1.2.0-gsl-2.3-ifort-16.0.2/include/fgsl/"
export FGSL_INC
FGSL_LIB=${LIB_PATH}"/FGSL/install/fgsl-1.2.0-gsl-2.3-ifort-16.0.2/lib/libfgsl.so"
export FGSL_LIB
export LD_LIBRARY_PATH=$FGSL_LIB:$LD_LIBRARY_PATH
FGSL_LIB_FLAG='-lfgsl'
export FGSL_LIB_FLAG
}


function GSL_for_COARSEAIR {
#GSL
GSL_DIR=${LIB_PATH}'/GSL/install/gsl-1.14-ifort-16.0.2/'
export GSL_DIR
export GSL_ROOT_DIR=$GSL_DIR
GSL_INC=${LIB_PATH}'/GSL/install/gsl-1.14-ifort-16.0.2/include/'
export GSL_INC
GSL_LIB=${LIB_PATH}'/GSL/install/gsl-1.14-ifort-16.0.2/lib/'
export GSL_LIB
GSL_LIB=${LIB_PATH}'/GSL/install/gsl-1.14-ifort-16.0.2/lib/'
export GSL_LIBS
export LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH
export PATH=${LIB_PATH}'/GSL/install/gsl-1.14-ifort-16.0.2/bin':$PATH
export GSL_FLAG='libgsl.so'
export GSL_CBLAS_FLAG='libgslcblas.so'
}

function GSL_for_SMUQ {
GSL_DIR=${LIB_PATH}'/GSL/install/gsl-2.3-ifort-16.0.2/'
export GSL_DIR
export GSL_ROOT_DIR=$GSL_DIR
GSL_INC=${LIB_PATH}'/GSL/install/gsl-2.3-ifort-16.0.2/include/gsl/'
export GSL_INC
GSL_LIB=${LIB_PATH}'/GSL/install/gsl-2.3-ifort-16.0.2/lib'
export GSL_LIB
gsl_LIBS=${LIB_PATH}'/GSL/install/gsl-2.3-ifort-16.0.2/lib'
export gsl_LIBS
export LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH
export PATH=${LIB_PATH}'/GSL/install/gsl-2.3-ifort-16.0.2/bin':$PATH
export GSL_FLAG='libgsl.so'
export GSL_CBLAS_FLAG='libgslcblas.so'
}
##################################################################################################################################


################################################################ CoarseAIR ##############################################################
function COARSEAIR_UPDATE {
#PLATO_LOAD
#KONIG_LOAD
INTEL_NEW_LOAD
GSL_for_COARSEAIR
cd $WORKSPACE_PATH/CoarseAIR/coarseair/
source ./scripts/building/BuildingInstalling.sh
}
###################################################################################################################################
