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

# Ex. of Call: bash RunKONIG.sh O3 10000 $WORKSPACE_PATH/Mars_Database/Run_0D/database/kinetics/ 1 1 1

System=${1}
echo "[RunMasterEq]: System                = "$System
TTran=${2}
echo "[RunMasterEq]: Transl. T             = "$TTran
PathToKinFldr=${3}
#PathToKinFldr=$WORKSPACE_PATH/Mars_Database/Run_0D/database/kinetics/
echo "[RunMasterEq]: Path To Kinetic Fldr? = "$PathToKinFldr
DissFlg=${4}
echo "[RunMasterEq]: Writing Dissociation? = "$DissFlg
InelFlg=${5}
echo "[RunMasterEq]: Writing Inelastic?    = "$InelFlg
ExchFlg=${6}
echo "[RunMasterEq]: How Many Exchanges?   = "$ExchFlg


echo 'Units=cm^3/s' > $PathToKinFldr/'KineticsTEMP'

iDiss=1
while [ ${iDiss} -le ${DissFlg} ]; do 
	echo "Adding Dissociation Kinetics to File "$PathToKinFldr/'KineticsTEMP'
	cat $PathToKinFldr/$System'Diss_'$TTran'K.dat' >> $PathToKinFldr/'KineticsTEMP'
	iDiss=$((iDiss+1))
done

iInel=1
while [ ${iInel} -le ${InelFlg} ]; do 
	echo "Adding Inelastic Kinetics to File "$PathToKinFldr/'KineticsTEMP'
	cat $PathToKinFldr/$System'Inel_'$TTran'K.dat' >> $PathToKinFldr/'KineticsTEMP'
	iInel=$((iInel+1))
done

iExch=1
while [ ${iExch} -le ${InelFlg} ]; do 
	echo "Adding Exchange Kinetics to File "$PathToKinFldr/'KineticsTEMP'
	cat $PathToKinFldr/$System'Exch_Type'$iExch'_'$TTran'K.dat' >> $PathToKinFldr/'KineticsTEMP'
	iExch=$((iExch+1))
done