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
    


def Write_Rates_Thermal(Syst, Temp, InputData):

    mkdirs( InputData.FinalFldr )    
    PathToFile = InputData.FinalFldr + '/' + InputData.SystNameLong + '_KTh.csv'
    print('    [Write_Rates_Thermal]: Writing Thermal Rates in File: ' + PathToFile )
    with open(PathToFile, 'w') as csvTermo:
        Line       = '# T,KDiss,KInel' 
        for iExch in range(2, Syst.NProcTypes):
            Line = Line + ',KExch' + str(iExch-1)
        Line = Line + '\n'
        csvTermo.write(Line)
        TempMat = np.concatenate( (np.expand_dims(np.array(Temp.TranVec, dtype='float'), axis=1), Syst.RatesTh), axis=1 )
        np.savetxt(csvTermo, TempMat, delimiter=',')
    csvTermo.close()

    
                
def Write_DissRates_Thermal(Syst, Temp, InputData):

    mkdirs( InputData.FinalFldr )    
    PathToFile = InputData.FinalFldr + '/' + Syst.Molecule[0].Name + '_KTh_Diss.csv'
    print('    [Write_DissRates_Thermal]: Writing Dissociation Thermal Rates in File: ' + PathToFile )
    with open(PathToFile, 'w') as csvTermo:
        Line    = '# T,KDiss\n' 
        csvTermo.write(Line)
        TempMat = np.concatenate( (np.expand_dims(np.array(Temp.TranVec, dtype='float'), axis=1), Syst.RatesTh[:,0]), axis=1 )
        np.savetxt(csvTermo, TempMat, delimiter=',')
    csvTermo.close()



def Write_PartFuncsAndEnergies(Syst, Temp, InputData):

    mkdirs( InputData.Kin.WriteFldr + '/thermo/' )    
    for iMol in range(Syst.NMolecules):

        for iT in Temp.iTVec:

            PathToFileOrig = InputData.Kin.WriteFldr + '/thermo/' + Syst.Molecule[iMol].Name + '_Format'
            PathToFile     = InputData.Kin.WriteFldr + '/thermo/' + Syst.Molecule[iMol].Name + '_' + str(Temp.TranVec[iT-1])
            DestTemp       = shutil.copyfile(PathToFileOrig, PathToFile)
            print('Copied File to: ', DestTemp)

            with open(PathToFile, 'a') as f:
                Line = 'NB_ENERGY_LEVELS = ' + str(Syst.Molecule[iMol].NBins) + '\n'
                f.write(Line)
                np.savetxt(f, np.transpose(np.array([Syst.Molecule[iMol].T[iT-1].Q, Syst.Molecule[iMol].T[iT-1].LevelEeV])), fmt='%.8e    %.8e')
            f.close()



def Write_Kinetics(Syst, Temp, InputData, iT):

    mkdirs( InputData.Kin.WriteFldr + '/kinetics/' )    

    if (InputData.Kin.WriteDiss_Flg):
        DissKinetics = InputData.Kin.WriteFldr + '/kinetics/' + Syst.Name + 'Diss_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
        csvkinetics  = open(DissKinetics, 'w')

        print('      [Write_Kinetics]: Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
        for iLevel in range(Syst.Molecule[0].NBins):

            if (Syst.T[iT-1].Proc[0].Rates[iLevel,0] > 0.0):
                ProcName = Syst.Molecule[0].Name + '(' + str(iLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
                Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % Syst.T[iT-1].Proc[0].Rates[iLevel,0]
                csvkinetics.write(Line)
        
        csvkinetics.close()


    if (InputData.Kin.WriteInel_Flg):
        InelKinetics = InputData.Kin.WriteFldr + '/kinetics/' + Syst.Name + 'Inel_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
        csvkinetics  = open(InelKinetics, 'w')

        print('      [Write_Kinetics]: Writing Inelastic: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
        InelFile     = InputData.Kin.ReadFldr  + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '.csv'

        for iLevel in range(Syst.Molecule[0].NBins):
            for jLevel in range(Syst.Molecule[0].NBins):

                if ((Syst.T[iT-1].Proc[1].Rates[iLevel,jLevel] > 0.0) and (Syst.Molecule[0].LevelEEh[iLevel] > Syst.Molecule[0].LevelEEh[jLevel]) ):
                    ProcName = Syst.Molecule[0].Name + '(' + str(iLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '(' + str(jLevel+1) + ')+' + Syst.Atom[2].Name
                    Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,5\n' % Syst.T[iT-1].Proc[1].Rates[iLevel,jLevel]
                    csvkinetics.write(Line)
                    
        csvkinetics.close()


    if (InputData.Kin.WriteExch_Flg):

        for iExch in range (2, Syst.NProcTypes):
            ExchKinetics = InputData.Kin.WriteFldr + '/kinetics/' + Syst.Name + 'Exch_Type' + str(iExch-1) + '_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
            csvkinetics  = open(ExchKinetics, 'w')

            print('      [Write_Kinetics]: Writing Exchange Nb. '+ str(iExch-1) + ': ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
            InelFile     = InputData.Kin.ReadFldr  + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '.csv'

            for iLevel in range(Syst.Molecule[0].NBins):
                for jLevel in range(Syst.Molecule[Syst.ExchtoMol[iExch-2]].NBins):

                    if ((Syst.T[iT-1].ProcExch[iExch-2].Rates[iLevel,jLevel] > 0.0) and (Syst.Molecule[0].LevelEEh[iLevel] > Syst.Molecule[Syst.ExchtoMol[iExch-2]].LevelEEh[jLevel]) ):
                        ProcName = Syst.Molecule[0].Name + '(' + str(iLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '(' + str(jLevel+1) + ')+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,6\n' % Syst.T[iT-1].ProcExch[iExch-2].Rates[iLevel,jLevel]
                        csvkinetics.write(Line)
                        
            csvkinetics.close()



def Write_Kinetics_FromOverall(Syst, Temp, InputData):

    mkdirs( InputData.Kin.WriteFldr + '/kinetics/' )    

    for iT in Temp.iTVec:
        print('\nTemperature Nb ', iT, '; T = ', Temp.TranVec[iT-1], 'K')

        if (InputData.Kin.WriteInel_Flg == True):
            InelKinetics = InputData.Kin.WriteFldr + '/kinetics/' + Syst.Name + 'Inel_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
            csvkinetics  = open(InelKinetics, 'w')

            print('  Writing Inelastic: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
            InelFile     = InputData.Kin.ReadFldr  + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '.csv'
            with open(InelFile) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                next(readCSV)
                for row in readCSV:

                    if (float(row[iT+1]) > 0.0):
                        ProcName = Syst.Molecule[0].Name + '(' + str(row[0]) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '(' + str(row[1]) + ')+' + Syst.Atom[2].Name
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,5\n' % float(row[iT+1])
                        csvkinetics.write(Line)
                
                csvfile.close()
                csvkinetics.close()


        if (InputData.Kin.WriteExch_Flg == True):

            for iExch in range (2, Syst.NProcTypes):
                print('  Writing Exchange: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name )
                ExchKinetics = InputData.Kin.WriteFldr + '/kinetics/' + Syst.Name + 'Exch_Type' + str(iExch-1) + '_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
                csvkinetics  = open(ExchKinetics, 'w')

                ExchFile     = InputData.Kin.ReadFldr + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name + '_Exch.csv'
                with open(ExchFile) as csvfile:
                    readCSV = csv.reader(csvfile, delimiter=',')
                    next(readCSV)
                    for row in readCSV:

                        if (float(row[iT+1]) > 0.0):
                            ProcName = Syst.Molecule[0].Name + '(' + str(row[0]) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '(' + str(row[1]) + ')+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name
                            Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,6\n' % float(row[iT+1])
                            csvkinetics.write(Line)
                
                    csvfile.close()
                    csvkinetics.close()


        if (InputData.Kin.WriteDiss_Flg == True):
            DissKinetics = InputData.Kin.WriteFldr + '/kinetics/' + Syst.Name + 'Diss_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
            csvkinetics  = open(DissKinetics, 'w')

            print('  Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
            DissFile = InputData.Kin.ReadFldr   + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '.csv'
            with open(DissFile) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                next(readCSV)
                for row in readCSV:
                    #print(row[:])

                    if (float(row[iT]) > 0.0):
                        ProcName = Syst.Molecule[0].Name + '(' + str(row[0]) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % float(row[iT])
                        csvkinetics.write(Line)
                
                csvfile.close()
                csvkinetics.close()


def Write_QSS(Syst, Temp, InputData, iT):

    PathToFile = InputData.FinalFldr + '/' + InputData.SystNameLong + '_KQSS.csv'
    print('    [Write]: Writing QSS Rates in File: ' + PathToFile )
    with open(PathToFile, 'a') as csvQSS:
        if (not path.exists(PathToFile) ):
            Line       = '# T,KDiss,KInel' 
            for iExch in range(2, Syst.NProcTypes):
                Line = Line + ',KExch' + str(iExch-1)
            Line = Line + '\n'
            csvQSS.write(Line)
        TempMat = np.transpose( np.expand_dims( np.concatenate( [np.array([Temp.TranVec[iT-1]], float), Syst.T[iT-1].QSS.Rate] ), axis=1 ) )
        np.savetxt(csvQSS, TempMat, delimiter=',')
    csvQSS.close()