export WORKSPACE_PATH='/home/venturi/WORKSPACE/'

###################################################### GSL & FGSL ################################################################
function GSL_OLD_LOAD {
GSL_DIR='/home/libo/LIBRARIES/gsl/gsl-1.14-gcc-7.3'
export GSL_DIR
GSL_INC=$GSL_DIR/'include'
export GSL_INC
GSL_LIB=$GSL_DIR/'lib'
export GSL_LIB
export LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH
GSL_INSTALL=$GSL_DIR
export GSL_INSTALL
export PATH=$GSL_INSTALL:$PATH
}
##################################################################################################################################


############################################################## CoarseAIR ############################################################
function COARSEAIR_UPDATE {
module purge
module load 3.2.10/moduleapps/gcc/7.3.0
module load 3.2.10/moduleapps/openmpi/3.0.0
export OMP_NUM_THREADS=32
module load 3.2.10/moduleslibs/openblas/0.2.20
module load 3.2.10/moduleslibs/sundials/3.0.0
#module load sundials-openmp/3.1.0
GSL_OLD_LOAD
cd $WORKSPACE_PATH/CoarseAIR/coarseair/
source ./scripts/building/BuildingInstalling.sh
}
##################################################################################################################################


########################################################## PLATO & KONIG #########################################################
function PLATO_SUBSTITUTE {
export PLATO_INSTALL_DIR=$WORKSPACE_PATH/'neqplasma_QCT/plato/install/plato-'$PLATO_DISTRIBUTION'-'$COMPILER'-'$COMPILER_VERSION
export PLATO_INC=$PLATO_INSTALL_DIR/'include'
export PLATO_LIB=$PLATO_INSTALL_DIR/'lib'
export PLATO_LIBS=$PLATO_LIB
export LD_LIBRARY_PATH=$PLATO_LIB:$LD_LIBRARY_PATH
export PLATO_EXECUTABLE=$PLATO_INSTALL_DIR/'bin/plato'
export PLATO_SOURCE_DIR=$WORKSPACE_PATH/'neqplasma_QCT/plato/plato-1.0'
export PLATO_DATABASE_DIR=$WORKSPACE_PATH/'neqplasma_QCT/database'
}

function PLATO_INSTALL {
cd $PLATO_SOURCE_DIR
make distclean
./autogen.sh
./configure --enable-optim --enable-omp --enable-parkin --enable-cache --prefix=${PLATO_INSTALL_DIR}
make all
make install
}


function PLATO_gnu_debug {
export PLATO_DISTRIBUTION='debug'
PLATO_SUBSTITUTE
}

function PLATO_gnu_release {
export PLATO_DISTRIBUTION='release'
PLATO_SUBSTITUTE
}

function PLATO_intel_debug {
export PLATO_DISTRIBUTION='debug'
PLATO_SUBSTITUTE
}

function PLATO_intel_release {
export PLATO_DISTRIBUTION='release'
PLATO_SUBSTITUTE
}

function PLATO_LOAD {
PLATO_gnu_release
}

function KONIG_SUBSTITUTE {
export KONIG_TYPE='orig'
export KONIG_INSTALL_DIR=$WORKSPACE_PATH/'neqplasma_QCT/konig/install/konig-'$KONIG_DISTRIBUTION'-'$KONIG_TYPE'-'$COMPILER'-'$COMPILER_VERSION
export KONIG_INC=$KONIG_INSTALL_DIR/'include'
export KONIG_LIB=$KONIG_INSTALL_DIR/'lib'
export KONIG_LIBS=$KONIG_LIB
export LD_LIBRARY_PATH=$KONIG_LIB:$LD_LIBRARY_PATH
export KONIG_EXECUTABLE=$KONIG_INSTALL_DIR/'bin/konig'
export KONIG_SOURCE_DIR=$WORKSPACE_PATH/'neqplasma_QCT/konig/konig'
export KONIG_DATABASE_DIR=$WORKSPACE_PATH/'neqplasma_QCT/database'
}

function KONIG_INSTALL {
cd $KONIG_SOURCE_DIR
make distclean
./autogen.sh
./configure --prefix=${KONIG_INSTALL_DIR} --enable-optim --enable-omp --with-oblas=${OBLAS_DIR} --with-plato=$PLATO_INSTALL_DIR
make all
make install
}

function KONIG_gnu_debug {
export KONIG_DISTRIBUTION='debug'
KONIG_SUBSTITUTE
}

function KONIG_gnu_release {
export KONIG_DISTRIBUTION='release'
KONIG_SUBSTITUTE
module load 3.2.10/moduleslibs/sundials/3.0.0
alias CVODE_BOX=$WORKSPACE_PATH/'neqplasma_QCT/cvode/cvode_examples/box/exec/box_'
alias CVODE_BOX_CG=$WORKSPACE_PATH/'neqplasma_QCT/cvode/cvode_examples/box_CG/exec/box_'
alias CVODE_BOX_2000K=$WORKSPACE_PATH/'neqplasma_QCT/cvode/cvode_examples/box_2000K/exec/box_'
alias CVODE_BOX_CG_2000K=$WORKSPACE_PATH/'neqplasma_QCT/cvode/cvode_examples/box_CG_2000K/exec/box_'
alias CVODE_BOX_3000K=$WORKSPACE_PATH/'neqplasma_QCT/cvode/cvode_examples/box_3000K/exec/box_'
alias CVODE_BOX_CG_3000K=$WORKSPACE_PATH/'neqplasma_QCT/cvode/cvode_examples/box_CG_3000K/exec/box_'
}

function KONIG_intel_debug {
export KONIG_DISTRIBUTION='debug'
KONIG_SUBSTITUTE
}

function KONIG_intel_release {
export KONIG_DISTRIBUTION='release'
KONIG_SUBSTITUTE
}

function KONIG_LOAD {
KONIG_gnu_release
}
################################################################################################################################## 
