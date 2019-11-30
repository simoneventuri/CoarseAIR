#!/bin/bash
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

Test1=1
Test2=1
Test3=1
Test4=1
Test5=1
Test6=1
Test7=1
Test8=1
Test9=1
Test10=0
Test11=0
Test12=1
Test13=1


export TEST_DIR=$(pwd)/../../test_New
mkdir -p  $TEST_DIR
export COARSEAIR_DTB_DIR=$(pwd)/../../dtb

scp -r ./TestN3 $TEST_DIR
scp ./SingleProcessor.sh $TEST_DIR/TestN3;  scp ./MultipleProcessors.sh $TEST_DIR/TestN3;  scp ./MultipleNodes.sh $TEST_DIR/TestN3


scp -r ./TestCO2 $TEST_DIR
scp ./SingleProcessor.sh $TEST_DIR/TestCO2; scp ./MultipleProcessors.sh $TEST_DIR/TestCO2; scp ./MultipleNodes.sh $TEST_DIR/TestCO2


### -------- TESTING N3 SYSTEM -------------------------------------------------------------------###
  echo; echo "  Generating Tests For N3 System ... "
  cd $TEST_DIR/TestN3
  
  export RUN_DIR=$(pwd)
  export COARSEAIR_WORKING_DIR=$(pwd)
  export COARSEAIR_INPUT_DIR=$(pwd)/input  
  
  
  ### ------ Testing StateToState ----------------------------------------------------------------###
    export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_StS
  
    #### --- Test 1 ------------------------------------------------------------------------------###
      if [ $Test1 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_StS_1.inp
        echo; echo "    Test 1 -> StS on 1 Processor @ 5000K, 10000K, 50000K for Level 1"
        ls
        bash SingleProcessor.sh #&>/dev/null
        echo "    Test 1 Generated! "
      fi
        
        
    #### --- Test 2 ------------------------------------------------------------------------------###
      if [ $Test2 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_StS_2.inp
        echo; echo "    Test 2 -> StS on 4 Processors @ 6000K, 20000K, 60000K for Level 5000"
        bash MultipleProcessors.sh #&>/dev/null
        echo "    Test 2 Generated! "
      fi
      
    #### --- Test 3 ------------------------------------------------------------------------------###
      if [ $Test3 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_StS_3.inp
        echo; echo "    Test 3 -> StS on 1 Processor @ 7000K, 30000K, 70000K for Level 9000"
        bash SingleProcessor.sh #&>/dev/null
        echo "    Test 3 Generated! "
      fi
            
  ### --------------------------------------------------------------------------------------------###   
    
    
  ### --- Testing VibrationalSpecific ------------------------------------------------------------###
    export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_VS
  
    #### --- Test 4 ------------------------------------------------------------------------------###
      if [ $Test4 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_VS_1.inp
        echo; echo "    Test 4 -> VibSpecific on 4 Processors @ 500K for vqn 1"
        bash MultipleProcessors.sh #&>/dev/null
        echo "    Test 4 Generated! "
      fi
      
    #### --- Test 5 ------------------------------------------------------------------------------###
      if [ $Test5 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_VS_2.inp
        echo; echo "    Test 5 -> VibSpecific on 1 Processor @ 2000K for vqn 80"
        bash SingleProcessor.sh #&>/dev/null
        echo "    Test 5 Generated! "
      fi
            
  ### --------------------------------------------------------------------------------------------###  
  
  
  ### --- Testing CoarseGrained ------------------------------------------------------------------###  
    export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_CG
    
    #### --- Test 6 ------------------------------------------------------------------------------###
      if [ $Test6 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_CG.inp
        echo; echo "    Test 6 -> CG60 on 4 Processors @ 5000K, 10000K, 20000K"
        bash MultipleProcessors.sh #&>/dev/null
        echo "    Test 6 Generated! "
      fi
      
  ### --------------------------------------------------------------------------------------------###  
  
  
###_______________________________________________________________________________________________###



### -------- TESTING CO2 SYSTEM ------------------------------------------------------------------###
  echo; echo "  Generating Tests For CO2 System ... "
  cd $TEST_DIR/TestCO2
  
  export RUN_DIR=$(pwd)
  export COARSEAIR_WORKING_DIR=$(pwd)
  export COARSEAIR_INPUT_DIR=$(pwd)/input  
  
  
  ### ------ Testing StateToState ----------------------------------------------------------------###
  
    #### --- Test 7 ------------------------------------------------------------------------------###
      if [ $Test7 -eq 1 ]; then
        export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_StS_1
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_StS_1.inp
        echo; echo "    Test 7 -> StS on 1 Processor @ 5000K, 10000K, 50000K for Level 7000, PES NASA_13A1"
        bash SingleProcessor.sh #&>/dev/null  
        echo "    Test 7 Generated! "
      fi
      
    #### --- Test 8 ------------------------------------------------------------------------------###
      if [ $Test8 -eq 1 ]; then
        export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_StS_2
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_StS_2.inp
        echo; echo "    Test 8 -> StS on 4 Processors @ 5000K, 10000K, 50000K for Level 7000, PES NASA_13A2"
        bash MultipleProcessors.sh #&>/dev/null
        echo "    Test 8 Generated! "
      fi
      
    #### --- Test 9 ------------------------------------------------------------------------------###
      if [ $Test9 -eq 1 ]; then
        export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_StS_3
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_StS_3.inp
        echo; echo "    Test 9 -> StS on 1 Processor @ 5000K, 10000K, 50000K for Level 7000, PES NASA_23A2"
        bash SingleProcessor.sh #&>/dev/null
        echo "    Test 9 Generated! "
      fi
            
  ### --------------------------------------------------------------------------------------------###   
    
    
  ### --- Testing VibrationalSpecific ------------------------------------------------------------###
    export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_VS
  
    #### --- Test 10 ------------------------------------------------------------------------------###
      if [ $Test10 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_VS_1.inp
        echo; echo "    Test 10 -> VibSpecific on 4 Processors @ 5000K, 10000K, 20000K"
        bash MultipleProcessors.sh #&>/dev/null  
        echo "    Test 10 Generated! "
      fi
      
    #### --- Test 11 -----------------------------------------------------------------------------###
      if [ $Test11 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_VS_2.inp
        echo; echo "    Test 11 -> VibSpecific on 1 Processor @ 5000K, 10000K, 20000K"
        bash SingleProcessor.sh #&>/dev/null        
        echo "    Test 11 Generated! "
      fi

  ### --------------------------------------------------------------------------------------------###  
  
  
  ### --- Testing CoarseGrained ------------------------------------------------------------------###  
    export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_CG
    
    #### --- Test 12 -----------------------------------------------------------------------------###
      if [ $Test12 -eq 1 ]; then
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_CG.inp
        echo; echo "    Test 12 -> CG60 on 4 Processors @ 5000K, 10000K, 20000K, NASA_13A1"
        bash MultipleProcessors.sh #&>/dev/null
        echo "    Test 12 Generated! "
      fi
      
  ### ---------------------------------------------------------------------------------------------###  
  
  
  ### --- Testing StateToState --------------------------------------------------------------------###  
  
    #### --- Test 13 ------------------------------------------------------------------------------###
      if [ $Test13 -eq 1 ]; then
        export COARSEAIR_OUTPUT_DIR=$(pwd)/BenchMark_StS_Merged
        export COARSEAIR_INPUT_FILE=$(pwd)/input/Test_StS_4.inp
        echo; echo "    Test 13 -> StS, Merging PESs @ 5000K, 10000K, 50000K for Level 7000"
        bash SingleProcessor.sh #&>/dev/null
        echo "    Test 13 Generated! "
      fi
  
  
###_______________________________________________________________________________________________###

echo; echo "REGRESSION TESTS GENERATED!!!"; echo

exit 0
