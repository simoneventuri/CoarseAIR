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


echo "Units=cm^3/s" > $PathToDtbFldr"/kinetics/KineticsTEMP_T"$TTran"K"

if [ $DissFlg -eq 1 ]; then
	echo "  [Initialize_0D_Database]: Adding Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
elif [ $DissFlg -eq 2 ]; then
	echo "  [Initialize_0D_Database]: Adding Corrected Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss_Corrected.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
elif [ $DissFlg -eq 3 ]; then
	echo "  [Initialize_0D_Database]: Adding VS Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss_VS.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"		
elif [ $DissFlg -eq 4 ]; then
	echo "  [Initialize_0D_Database]: Adding Phys-Based Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss_Phys_"${NBins}"Bins.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"	
elif [ $DissFlg -eq 5 ]; then
	echo "  [Initialize_0D_Database]: Adding Fitted Phys-Based Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss_Phys_Fitted_"${NBins}"Bins.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"	
elif [ $DissFlg -eq 6 ]; then
	echo "  [Initialize_0D_Database]: Adding Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	echo "  [Initialize_0D_Database]: MANINDER'S CASE: NO CORRECTION and NO RECOMBINATION"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
elif [ $DissFlg -eq 7 ]; then
	echo "  [Initialize_0D_Database]: Adding Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	echo "  [Initialize_0D_Database]: MANINDER'S CASE: NO CORRECTION but RECOMBINATION"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"	
elif [ $DissFlg -eq 11 ]; then
	echo "  [Initialize_0D_Database]: Adding Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	echo "  [Initialize_0D_Database]: NO CORRECTION and NO RECOMBINATION"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"	
elif [ $DissFlg -eq 12 ]; then
	echo "  [Initialize_0D_Database]: Adding Corrected Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	echo "  [Initialize_0D_Database]: NO RECOMBINATION"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Diss.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"	
fi

# iDissInel=${DissInelFlg}
# if [ ${iDissInel} -eq 1 ]; then
# 	echo "  [Initialize_0D_Database]: Adding Inelastic + Dissociation Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
# 	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/DissInel.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
# # elif [ ${iDissInel} -eq 2 ]; then
# # 	echo "  [Initialize_0D_Database]: Adding Window-Averaged Inelastic Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
# # 	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Inel_WindAvrg.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
# fi

iInel=${InelFlg}
if [ ${iInel} -eq 1 ]; then
	echo "  [Initialize_0D_Database]: Adding Inelastic Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Inel.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
elif [ ${iInel} -eq 2 ]; then
	echo "  [Initialize_0D_Database]: Adding Window-Averaged Inelastic Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Inel_WindAvrg.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
fi

iExch=${ExchFlg1}
if [ ${iExch} -eq 1 ]; then
	echo "  [Initialize_0D_Database]: Adding Exchange Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Exch_Type1.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
elif [ ${iExch} -eq 2 ]; then
	echo "  [Initialize_0D_Database]: Adding Window-Averaged Exchange Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Exch_Type1_WindAvrg.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
fi

iExch=${ExchFlg2}
if [ ${iExch} -eq 1 ]; then
	echo "  [Initialize_0D_Database]: Adding Exchange Kinetics to File "$PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
	cat $PathToDtbFldr"/kinetics/"${System}${FldrName}"/T"$TTran"K/Exch_Type2.dat" >> $PathToDtbFldr/"/kinetics/KineticsTEMP_T"$TTran"K"
fi

for Molecule in "${Molecule_vec[@]}"; do :

	echo "  [Initialize_0D_Database]: Copying Thermo File "$PathToDtbFldr/"/thermo/"${Molecule}"_T"${TTran}"K"
	cat $PathToDtbFldr/"/thermo/"${System}"/"${Molecule}"_T"${TTran}"K" > $PathToDtbFldr/"/thermo/"${Molecule}"_T"${TTran}"K"
#	cat $PathToDtbFldr/"/thermo/"${System}${FldrName}"/"${Molecule}"_T"${TTran}"K" > $PathToDtbFldr/"/thermo/"${Molecule}"_T"${TTran}"K"

	echo "  [Initialize_0D_Database]: Copying Initial Mole Fraction File "$PathToDtbFldr/"/thermo/"${Molecule}"_T"${T0}"K"
	cat $PathToDtbFldr/"/thermo/"${System}"/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat" > $PathToDtbFldr/"/thermo/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat"
#	cat $PathToDtbFldr/"/thermo/"${System}${FldrName}"/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat" > $PathToDtbFldr/"/thermo/"${Molecule}"_InitialMoleFracs_T"${T0}"K.dat"

done

}