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
function Initialize_0D_Database {
# Ex. of Call: bash Initialize_0D_Database.sh O3 10000 $WORKSPACE_PATH/Mars_Database/Run_0D/database/kinetics/ 1 1 1

#System=${1}
#echo "  [Initialize_0D_Database]: System                = "$System
#TTran=${2}
#echo "  [Initialize_0D_Database]: Transl. T             = "$TTran
#PathToDtbFldr=${3}
#PathToDtbFldr=$WORKSPACE_PATH/Mars_Database/Run_0D/database/kinetics/
#echo "  [Initialize_0D_Database]: Path To Kinetic Fldr? = "$PathToDtbFldr
#DissFlg=${4}
#echo "  [Initialize_0D_Database]: Writing Dissociation? = "$DissFlg
#InelFlg=${5}
#echo "  [Initialize_0D_Database]: Writing Inelastic?    = "$InelFlg
#ExchFlg=${6}
#echo "  [Initialize_0D_Database]: How Many Exchanges?   = "$ExchFlg


echo "Units=cm^3/s" > $PathToDtbFldr"/kinetics/KineticsTEMP"

iDiss=1
while [ ${iDiss} -le ${DissFlg} ]; do 
	echo "  [Initialize_0D_Database]: Adding Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP"
	cat $PathToDtbFldr"/kinetics/"$System"Diss_"$TTran"K.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP"
	iDiss=$((iDiss+1))
done

iInel=1
while [ ${iInel} -le ${InelFlg} ]; do 
	echo "  [Initialize_0D_Database]: Adding Inelastic Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP"
	cat $PathToDtbFldr"/kinetics/"$System"Inel_"$TTran"K.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP"
	iInel=$((iInel+1))
done

iExch=${ExchFlg1}
while [ ${iExch} -gt 0 ] || [ ${iExch} -le ${ExchFlg2} ]; do 
	echo "  [Initialize_0D_Database]: Adding Exchange Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP"
	cat $PathToDtbFldr"/kinetics/"$System"Exch_Type"$iExch"_"$TTran"K.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP"
	iExch=$((iExch+1))
done

}