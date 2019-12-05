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

# ------------------------------------------------------------------------------------------------------------- DeriveQuantities #
function DeriveQuantities {

  echo "  [DeriveQuantities]: COARSEAIR_SOURCE_DIR  = "${COARSEAIR_SOURCE_DIR}
  echo "  [DeriveQuantities]: COARSEAIR_OUTPUT_DIR  = "${COARSEAIR_OUTPUT_DIR}
  echo "  [DeriveQuantities]: KONIG_SOURCE_DIR      = "${KONIG_SOURCE_DIR}
  echo "  [DeriveQuantities]: KONIG_DATABASE_DIR    = "${KONIG_DATABASE_DIR}
  echo "  [DeriveQuantities]: NLevels1              = "${NLevels1}
  echo "  [DeriveQuantities]: MinLevel1             = "${MinLevel1}
  echo "  [DeriveQuantities]: MaxLevel1             = "${MaxLevel1}
  echo "  [DeriveQuantities]: NLevels2              = "${NLevels2}
  echo "  [DeriveQuantities]: MinLevel2             = "${MinLevel2}
  echo "  [DeriveQuantities]: MaxLevel2             = "${MaxLevel2}
  echo "  [DeriveQuantities]: StochPESFlg           = "${StochPESFlg}
  echo "  [DeriveQuantities]: NPESs                 = "${NPESs}
  echo "  [DeriveQuantities]: iPESStart             = "${iPESStart}
  echo "  [DeriveQuantities]: Molecule1             = "${Molecule1}
  echo "  [DeriveQuantities]: Tran                  = "${Tran}
  echo "  [DeriveQuantities]: System                = "${System}
  echo "  [DeriveQuantities]: RunExtCodeFlg         = "${RunExtCodeFlg}

  DeriveQuantitiesCommand="coarseair-derivequantities.x "${SilencingString}

  if [[ ${RunExtCodeFlg} -gt 2 ]]; then
    echo "  [DeriveQuantities]: Creating Directories for HEGEL"

    for Tran in "${Tran_vec[@]}"; do :
      echo '[DeriveQuantities] '
      if [ ${TranFlg} -eq 0 ]; then 
        echo "[DeriveQuantities]: ===   Translational Energy           = " ${Tran} " ============================ "
        #for Tint in "${Tint_Vec[@]}"; do :
        Tint=0.0
      else
        echo "[DeriveQuantities]: ===   Translational Temperature      = " ${Tran} " ============================ "
        #for Tint in "${Tint_Vec[@]}"; do :
        Tint=${Tran}
      fi
      echo "[DeriveQuantities]: ===== Internal Temperature           = " ${Tint} " ======================== "
      echo '[DeriveQuantities] '   
    done 

    
    if [ ${StochPESFlg} -eq 0 ]; then 
      NPESs=1
      iPES=1
    else
      iPES=${iPESStart}
    fi
    while [ ${iPES} -le ${NPESs} ]
    do
      echo "  [DeriveQuantities]: Computing Arrhenius Fittings for PES "${iPES}". Command: "${DeriveQuantitiesCommand}

      cd ${COARSEAIR_OUTPUT_DIR}
      mkdir -p RunHegel
      #scp   -r ${HEGEL_SOURCE_DIR}/input       ./RunHegel/
      mkdir -p ./RunHegel/database
      mkdir -p ./RunHegel/database/kinetics
      scp   -r ${HEGEL_DATABASE_DIR}/mixture   ./RunHegel/database
      scp   -r ${HEGEL_DATABASE_DIR}/thermo    ./RunHegel/database
      scp   -r ${HEGEL_DATABASE_DIR}/transfer  ./RunHegel/database
      scp   -r ${HEGEL_DATABASE_DIR}/transport ./RunHegel/database
      mkdir -p ./RunHegel/postprocessing

      if [ ${StochPESFlg} -eq 1 ]; then 
        mv ${COARSEAIR_OUTPUT_DIR}/RunKonig ${COARSEAIR_OUTPUT_DIR}/RunKonig_${iPES}
        cd ${COARSEAIR_OUTPUT_DIR}/RunHegel_${iPES}
      else
        cd ${COARSEAIR_OUTPUT_DIR}/RunHegel
      fi
      
      eval ${DeriveQuantitiesCommand} 0 ${iPES} ${NLevels1} ${MinLevel1} ${MaxLevel1} ${NLevels2} ${MinLevel2} ${MaxLevel2}

      iPES=$((iPES+1))
    done
    
    #scp -p ${COARSEAIR_SOURCE_DIR}/extra/PostHEGEL/Exec/PostHEGEL.m ./postprocessing

  elif [[ ${RunExtCodeFlg} -gt 0 ]]; then

    echo "  [DeriveQuantities]: Creating Directories for KONIG"
    cd ${COARSEAIR_OUTPUT_DIR}
    mkdir -p RunKonig
    scp   -r ${KONIG_SOURCE_DIR}/input       ./RunKonig/
    mkdir -p ./RunKonig/database
    mkdir -p ./RunKonig/database/kinetics
    scp   -r ${KONIG_DATABASE_DIR}/mixture   ./RunKonig/database
    scp   -r ${KONIG_DATABASE_DIR}/thermo    ./RunKonig/database
    scp   -r ${KONIG_DATABASE_DIR}/transfer  ./RunKonig/database
    scp   -r ${KONIG_DATABASE_DIR}/transport ./RunKonig/database
    mkdir -p ./RunKonig/postprocessing

    
    if [ ${StochPESFlg} -eq 0 ]; then 
      NPESs=1
      iPES=1
    else
      iPES=${iPESStart}
    fi
    
    while [ ${iPES} -le ${NPESs} ]
    do
      echo "  [DeriveQuantities]: Running 0D Simulation for PES "${iPES}". Command: "${DeriveQuantitiesCommand}
      
      cd ${COARSEAIR_OUTPUT_DIR}/RunKonig
      eval ${DeriveQuantitiesCommand} 0 ${iPES} ${NLevels1} ${MinLevel1} ${MaxLevel1} ${NLevels2} ${MinLevel2} ${MaxLevel2}
      
      if [[ ${RunExtCodeFlg} -eq 2 ]]; then
        if [ ${StochPESFlg} -eq 1 ]; then 
          rm -rf ${COARSEAIR_OUTPUT_DIR}/RunKonig/database/kinetics/${System}
          mv ${COARSEAIR_OUTPUT_DIR}/RunKonig/output/T_${Tran%.*}/output/box.dat              ${COARSEAIR_OUTPUT_DIR}/RunKonig/output/T_${Tran%.*}/output/box.dat.${iPES}
          mv ${COARSEAIR_OUTPUT_DIR}/RunKonig/output/T_${Tran%.*}/output/pop_${Molecule1}.dat ${COARSEAIR_OUTPUT_DIR}/RunKonig/output/T_${Tran%.*}/output/pop_${Molecule1}.dat.${iPES}
        fi
      else
        if [ ${StochPESFlg} -eq 1 ]; then 
          mv ${COARSEAIR_OUTPUT_DIR}/RunKonig/database/kinetics/${System} ${COARSEAIR_OUTPUT_DIR}/RunKonig/database/kinetics/${System}.${iPES}
        fi
      fi
      
      iPES=$((iPES+1))
    done
    if [[ ${RunExtCodeFlg} -eq 2 ]]; then
      scp -p ${COARSEAIR_SOURCE_DIR}/extra/PostKONIG/Exec/PostKONIG.m ./postprocessing
    fi

  else

    echo "  [DeriveQuantities]: Creating Directories for running DeriveQuantities"
    cd ${COARSEAIR_OUTPUT_DIR}
    mkdir -p DeriveQuantities
    cd ${COARSEAIR_OUTPUT_DIR}/DeriveQuantities
    
    eval ${DeriveQuantitiesCommand} 0 1 ${NLevels1} ${MinLevel1} ${MaxLevel1} ${NLevels2} ${MinLevel2} ${MaxLevel2}

  fi
  echo " "
  echo "  [DeriveQuantities]: -> Done with DeriveQuantities"
  echo "  [DeriveQuantities]: ------------------------------------------------------------ "
  echo " "
}
#================================================================================================================================#
