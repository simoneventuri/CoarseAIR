#===============================================================================================================
# 
# Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
# 
# Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
#
# Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of the 
# Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU Lesser General Public License for more details. 
# 
# You should have received a copy of the GNU Lesser General Public License along with this library; 
# if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
# 
#---------------------------------------------------------------------------------------------------------------
#===============================================================================================================

function COARSEAIR_SUBSTITUTE {
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     Machine=Linux;;
    Darwin*)    Machine=Mac;;
    CYGWIN*)    Machine=Cygwin;;
    MINGW*)     Machine=MinGw;;
    *)          Machine="UNKNOWN:${unameOut}"
esac
echo "Machine         : "$Machine

export WORKSPACE_PATH=$WORKSPACE_PATH #$(pwd)/../../
echo "Path to WORKSPACE: "$WORKSPACE_PATH
if [ ${FC} == "gfortran" ]; then
  export COMPILER_CGQCT="gnu"
  export COMPILER_CGQCT_VER=$(gcc -dumpversion)
  export USE_OPENBLAS_FLAG="YES"
  export NB_PROC=8
elif [ ${FC} == "ifort" ]; then
  export COMPILER_CGQCT="intel"
  export COMPILER_CGQCT_VER=18.0.3 #$COMPILER_VERSION_MICRO      #### PLEIADES
  #export COMPILER_CGQCT_VER=19.1.1 #$COMPILER_VERSION_MICRO     #### LAB
  export USE_OPENBLAS_FLAG="NO"
  export NB_PROC=8
fi
echo "Compiler        : "$COMPILER_CGQCT
echo "Compiler Version: "$COMPILER_CGQCT_VER
export COARSEAIR_LIBRARY_NAME="coarseair"
export COARSEAIR_VERSION="1.1"
export COARSEAIR_BUILD_CONFIG=$COARSEAIR_LIBRARY_NAME"-"$COARSEAIR_VERSION"-"$COARSEAIR_DISTRIBUTION"-"$COMPILER_CGQCT"-"$COMPILER_CGQCT_VER
export COARSEAIR_SOURCE_DIR=$WORKSPACE_PATH/"CoarseAIR/coarseair"
export COARSEAIR_BUILD_DIR=$WORKSPACE_PATH/"CoarseAIR/build"/$COARSEAIR_BUILD_CONFIG
export COARSEAIR_INSTALL_DIR=$WORKSPACE_PATH/"CoarseAIR/install"/$COARSEAIR_BUILD_CONFIG
export COARSEAIR_LIBRARY_DIR=$COARSEAIR_INSTALL_DIR/"lib"
export COARSEAIR_MODULES_DIR=$COARSEAIR_INSTALL_DIR/"mod"
export COARSEAIR_CMAKE_DIR=$COARSEAIR_INSTALL_DIR/"cmake"
export COARSEAIR_EXECUTABLE_DIR=$COARSEAIR_INSTALL_DIR/"bin"
export COARSEAIR_DATABASE_DIR=$COARSEAIR_SOURCE_DIR/"DATABASE"
export COARSEAIR_RUNS_DIR=$WORKSPACE_PATH/"CoarseAIR/runs"
export COARSEAIR_INPUT_DIR=$COARSEAIR_RUNS_DIR/"input"
export COARSEAIR_INPUT_FILE=$COARSEAIR_INPUT_DIR/"CoarseAIR.inp"
export COARSEAIR_OUTPUT_DIR=$COARSEAIR_RUNS_DIR/"Test"
export COARSEAIR_BIN_OUTPUT_DIR=$COARSEAIR_OUTPUT_DIR
export PATH=$COARSEAIR_EXECUTABLE_DIR:$PATH
export UseGSL='YES'
}

function COARSEAIR_INSTALL {
mkdir -p $COARSEAIR_BUILD_DIR
cd $COARSEAIR_BUILD_DIR
if [ ${Machine} == "Mac" ]; then
  cmake ../../coarseair -DCMAKE_BUILD_TYPE=$COARSEAIR_DISTRIBUTION -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_INSTALL_PREFIX=../../install -DGSL_INCLUDE_DIR=$GSL_INC -DGSL_LIBRARIES=$GSL_LIB -DGSL_ROOT_DIR=$GSL_DIR -DUSE_OPENBLAS=$USE_OPENBLAS_FLAG -DFGSL_INC=$FGSL_INC -DFGSL_LIB=$FGSL_LIB -DOpenBLAS_LIB_DIR=$OBLAS_LIB -DOpenBLAS_INCLUDE_DIR=$OBLAS_INC -DHDF5_LIBRARIES=$HDF5_LIB/libhdf5_fortran.dylib -DHDF5_INCLUDE_DIRS=$HDF5_INC -DHDF5_ROOT=$HDF5_DIR -DHDF5_Fortran_INCLUDE_DIRS=$HDF5_INC
else
  cmake ../../coarseair -DCMAKE_BUILD_TYPE=$COARSEAIR_DISTRIBUTION -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_INSTALL_PREFIX=../../install -DGSL_INCLUDE_DIR=$GSL_INC -DGSL_LIBRARIES=$GSL_LIB/libgsl.so -DGSL_ROOT_DIR=$GSL_DIR -DUSE_OPENBLAS=$USE_OPENBLAS_FLAG -DFGSL_INC=$FGSL_INC -DFGSL_LIB=$FGSL_LIB -DOpenBLAS_LIB_DIR=$OBLAS_LIB -DOpenBLAS_INCLUDE_DIR=$OBLAS_INC -DHDF5_LIBRARIES=$HDF5_LIB/libhdf5_fortran.dylib -DHDF5_INCLUDE_DIRS=$HDF5_INC -DHDF5_ROOT=$HDF5_DIR -DHDF5_Fortran_INCLUDE_DIRS=$HDF5_INC
fi
make -j $NB_PROC all install
cd $COARSEAIR_SOURCE_DIR 
}

function COARSEAIR_TEST {
cd $COARSEAIR_SOURCE_DIR/"scripts/test"
bash RegressionTest.sh
cd $COARSEAIR_SOURCE_DIR 
}

function COARSEAIR_debug {
export COARSEAIR_DISTRIBUTION='debug'
COARSEAIR_SUBSTITUTE
}


function COARSEAIR_release {
export COARSEAIR_DISTRIBUTION='release'
COARSEAIR_SUBSTITUTE
}
