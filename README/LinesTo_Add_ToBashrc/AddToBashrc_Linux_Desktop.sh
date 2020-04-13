export WORKSPACE_PATH='/home/venturi/WORKSPACE/'


########################################################## COMPILERS #############################################################
#GNU
function LOAD_GNU {
export COMPILER_VER=$GCC_VER
export FORTRAN_COMPILER=gfortran  
export FC=gfortran
export F77=gfortran
export CMAKE_Fortran_COMPILER=gfortran
export CC=gcc
export CXX=g++
}

function GNU_93_LOAD {
export GCC_DIR='/home/venturi/APPLICATIONS/GCC/gcc-9.3.0-install/'
export GCC_INC=$GCC_DIR/'include'
export GCC_LIB=$GCC_DIR/'lib64/'
export LD_LIBRARY_PATH=$GCC_LIB:$LD_LIBRARY_PATH
export GCC_INSTALL=$GCC_DIR/'bin'
export PATH=$GCC_INSTALL:$PATH
export GCC_VER='gnu-9.3.0'
LOAD_GNU
alias gcc=$GCC_INSTALL/'gcc'
alias cc=$GCC_INSTALL/'cc'
alias g++=$GCC_INSTALL/'g++'
alias c++=$GCC_INSTALL/'c++'
}
##################################################################################################################################


########################################################## LIBRARIES #############################################################
# OpenBLAS
OBLAS_DIR='/home/venturi/LIBRARIES/OPENBLAS/OpenBLAS-0.3.9-gcc-9.3.0/opt/OpenBLAS/'
export OBLAS_DIR
OBLAS_INC=$OBLAS_DIR/'include'
export OBLAS_INC
OBLAS_LIB=$OBLAS_DIR/'lib'
export OBLAS_LIB
export LD_LIBRARY_PATH=$OBLAS_LIB:$LD_LIBRARY_PATH
export OPENBLAS_NUM_THREADS=4
#################################################################################################################################


############################################################ CoarseAIR ############################################################
function COARSEAIR_UPDATE {
GNU_93_LOAD
cd $WORKSPACE_PATH/CoarseAIR/coarseair/
source ./scripts/building/BuildingInstalling.sh
}
##################################################################################################################################