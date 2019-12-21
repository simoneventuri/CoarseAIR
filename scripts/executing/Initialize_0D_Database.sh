#!/bin/bash

System=${1}
TTran=${2}
DissFlg=${3}
InelFlg=${4}
ExchFlg=${5}

PathToKinFldr=$WORKSPACE_PATH/Mars_Database/Run_0D/database/kinetics/

echo 'Units=cm^3/s' > $PathToKinFldr/KineticsTEMP

iDiss=1
while [ ${iDiss} -le ${DissFlg} ]; do 
	cat $PathToKinFldr/$System'Diss_'$TTran'K.dat' >> $PathToKinFldr/KineticsTEMP
	iDiss=$((iDiss+1))
done

iInel=1
while [ ${iInel} -le ${InelFlg} ]; do 
	cat $PathToKinFldr/$System'Inel_'$TTran'K.dat' >> $PathToKinFldr/KineticsTEMP\
	iInel=$((iInel+1))
done

iExch=1
while [ ${iInel} -le ${InelFlg} ]; do 
	cat $PathToKinFldr/$System'Exch_'$TTran'K.dat' >> $PathToKinFldr/KineticsTEMP
	iExch=$((iExch+1))
done