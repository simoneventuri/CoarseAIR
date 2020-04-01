##==============================================================================================================
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
##==============================================================================================================
import numpy as np
import pandas
import csv
import shutil
import os, errno
import os.path
from os import path

def mkdirs(newdir, mode=0o777):
    os.makedirs(newdir, mode, exist_ok=True)
    



def Write_QSS(Syst, Temp, InputData, iT):

    PathToFile = InputData.FinalFldr + '/KQSS.csv'
    print('    [Write_QSS]: Writing QSS Rates in File: ' + PathToFile )
    if (not path.exists(PathToFile) ):
        WriteFlg = True
    else:
        WriteFlg = False        
    with open(PathToFile, 'a') as csvQSS:
        if (WriteFlg):
            Line       = '# T,KDiss,KInel' 
            for iExch in range(2, Syst.NProcTypes):
                Line = Line + ',KExch' + str(iExch-1)
            Line = Line + ',tIni,tFin,ActualQSS?\n'
            csvQSS.write(Line)
        RealFlg = 1.0
        if (Syst.T[iT-1].QSS.Rate[0] >= Syst.RatesTh[0,iT-1]):
            RealFlg = 0.0
        TempVec = np.array([Syst.T[iT-1].QSS.Time[0], Syst.T[iT-1].QSS.Time[1], RealFlg], float)
        TempMat = np.transpose( np.expand_dims( np.concatenate( [np.array([Temp.TranVec[iT-1]], float), Syst.T[iT-1].QSS.Rate, TempVec] ), axis=1 ) )
        np.savetxt(csvQSS, TempMat, delimiter=',')
    csvQSS.close()



def Write_Arrhenius_Diss(Syst, InputData, Ixd, Coeffs, MaxEntOrPlato, csvkinetics):
    
    if (csvkinetics == 0):
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName )
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName + '/Arrhenius/' )    
        TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName + '/Arrhenius/'
        if (MaxEntOrPlato == 1):
            DissKinetics = TempFldr + '/Dh.dat' 
            print('      [Write_Arrhenius_Diss]: Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
            csvkinetics  = open(DissKinetics, 'w')
            Line     = '#LevelId,A1,A2,A3\n'
            csvkinetics.write(Line)    
        elif (MaxEntOrPlato == 2):
            if (InputData.Kin.CorrFactor != 1.0):
                DissKinetics = TempFldr + '/Diss_Corrected.dat' 
                print('      [Write_Arrhenius_Diss]: Writing Corrected Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
            else:
                DissKinetics = TempFldr + '/Diss.dat' 
                print('      [Write_Arrhenius_Diss]: Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
            csvkinetics  = open(DissKinetics, 'w')
    else:
        if (MaxEntOrPlato == 1):
            Line     = '%i,%.10e,%.10e,%.10e\n' % ((Ixd+1), float(Coeffs[0]), float(Coeffs[1]), float(Coeffs[2])) 
            csvkinetics.write(Line)       
        else:
            ProcName = Syst.Molecule[0].Name + '(' + str(Ixd+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
            Line     = ProcName + ':%.4e,%.4e,%.4e,2\n' % (float(Coeffs[0]), float(Coeffs[1]), float(Coeffs[2]))
            csvkinetics.write(Line)

    return csvkinetics



def Write_Arrhenius_Inel(Syst, InputData, Ixd, Coeffs, MaxEntOrPlato, iExch, csvkinetics):
    
    if (csvkinetics == 0):
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName )
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName + '/Arrhenius/' )    
        TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName + '/Arrhenius/'
        if (MaxEntOrPlato == 1):
            if (iExch == 0):
                InelKinetics = TempFldr + '/EXh_WOExch.dat' 
                print('      [Write_Arrhenius_Inel]: Writing Inelastic Processes, WITHOUT Exchange: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
            elif (iExch == 1):
                InelKinetics = TempFldr + '/EXh.dat' 
                print('      [Write_Arrhenius_Inel]: Writing Inelastic Processes, WITH Exchange: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
            csvkinetics  = open(InelKinetics, 'w')
            Line     = '#InLevelId,FinLevelId,A1,A2,A3\n'
            csvkinetics.write(Line)    
        elif (MaxEntOrPlato == 2):
            InelKinetics = TempFldr + '/Inel.dat' 
            print('      [Write_Arrhenius_Inel]: Writing Inelastic Processes: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
            csvkinetics  = open(InelKinetics, 'w')
    else:
        if (MaxEntOrPlato == 1):
            Line     = '%i,%i,%.10e,%.10e,%.10e\n' % ((Ixd[0]+1), (Ixd[1]+1), float(Coeffs[0]), float(Coeffs[1]), float(Coeffs[2])) 
            csvkinetics.write(Line)       
        else:
            ProcName = Syst.Molecule[0].Name + '(' + str(Ixd[0]+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '(' + str(rIxd[1]+1) + ')+' + Syst.Atom[2].Name
            Line     = ProcName + ':%.4e,%.4e,%.4e,5\n' % (float(Coeffs[0]), float(Coeffs[1]), float(Coeffs[2]))
            csvkinetics.write(Line)

    return csvkinetics



def Write_Arrhenius_Exch(Syst, InputData, Ixd, Coeffs, iExch, csvkinetics):
    
    if (csvkinetics == 0):
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName )
        mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName + '/Arrhenius/' )    
        TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + InputData.Kin.Groups.FldrName + '/Arrhenius/'
        ExchKinetics = TempFldr + 'Exch_Type' + str(iExch-1) + '.dat'  
        print('      [Write_Arrhenius_Inel]: Writing Exchange Processes Nb. '+ str(iExch-1) + ': ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name  + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name  )
        csvkinetics  = open(ExchKinetics, 'w')
    else:    
        ProcName = Syst.Molecule[0].Name + '(' + str(iLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '(' + str(jLevel+1) + ')+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name
        Line     = ProcName + ':%.4e,%.4e,%.4e,6\n' % (float(Coeffs[0]), float(Coeffs[1]), float(Coeffs[2]))
        csvkinetics.write(Line)

    return csvkinetics