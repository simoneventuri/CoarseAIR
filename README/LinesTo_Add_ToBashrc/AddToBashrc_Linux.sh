export WORKSPACE_PATH='/home/venturi/WORKSPACE/'



############################################################## CoarseAIR ############################################################
function COARSEAIR_UPDATE {
module purge
module load 3.2.10/moduleapps/gcc/7.3.0
module load 3.2.10/moduleapps/openmpi/3.0.0
export OMP_NUM_THREADS=32
module load 3.2.10/moduleslibs/openblas/0.2.20
#module load 3.2.10/moduleslibs/sundials/3.0.0
#module load sundials-openmp/3.1.0
#GSL_OLD_LOAD
cd $WORKSPACE_PATH/CoarseAIR/coarseair/
source ./scripts/building/BuildingInstalling.sh
}
##################################################################################################################################