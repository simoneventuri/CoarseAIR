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
    PathToFile = InputData.FinalFldr + '/KTh.csv'
    print('    [Write_Rates_Thermal]: Writing Thermal Rates in File: ' + PathToFile )
    if (not path.exists(PathToFile) ):
        WriteFlg = True
    else:
        WriteFlg = False        
    with open(PathToFile, 'a') as csvTermo:
        if (WriteFlg):
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
    PathToFile = InputData.FinalFldr + '/KTh_Diss.csv'
    print('    [Write_DissRates_Thermal]: Writing Dissociation Thermal Rates in File: ' + PathToFile )
    if (not path.exists(PathToFile) ):
        WriteFlg = True
    else:
        WriteFlg = False        
    with open(PathToFile, 'a') as csvTermo:
        if (WriteFlg):
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
            PathToFile     = InputData.Kin.WriteFldr + '/thermo/' + Syst.Molecule[iMol].Name + '_' + Syst.NameLong + '_' + str(Temp.TranVec[iT-1])
            DestTemp       = shutil.copyfile(PathToFileOrig, PathToFile)
            print('Copied File to: ', DestTemp)

            with open(PathToFile, 'a') as f:
                Line = 'NB_ENERGY_LEVELS = ' + str(Syst.Molecule[iMol].NBins) + '\n'
                f.write(Line)
                np.savetxt(f, np.transpose(np.array([Syst.Molecule[iMol].T[iT-1].Q, Syst.Molecule[iMol].T[iT-1].LevelEeV])), fmt='%.8e    %.8e')
            f.close()



def Write_Kinetics(Syst, Temp, InputData, iT):

    mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
    mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong )
    mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + '/T' + str(Temp.TranVec[iT-1]) + 'K/' )    
    TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + '/T' + str(Temp.TranVec[iT-1]) + 'K/'

    if (InputData.Kin.WriteDiss_Flg):
        if (InputData.Kin.CorrFactor != 1.0):
            DissKinetics = TempFldr + '/Diss_Corrected.dat' 
            print('      [Write_Kinetics]: Writing Corrected Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
        else:
            DissKinetics = TempFldr + '/Diss.dat' 
            print('      [Write_Kinetics]: Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
        csvkinetics  = open(DissKinetics, 'w')

        for iLevel in range(Syst.Molecule[0].NBins):

            if (Syst.T[iT-1].Proc[0].Rates[iLevel,0] > 0.0):
                ProcName = Syst.Molecule[0].Name + '(' + str(iLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
                Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % float(Syst.T[iT-1].Proc[0].Rates[iLevel,0])
                csvkinetics.write(Line)
        
        csvkinetics.close()


    if (InputData.Kin.WriteInel_Flg):
        InelKinetics = TempFldr + '/Inel.dat' 
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
            print('      [Write_Kinetics]: iExch =  ' + str(iExch-1) )

            ExchKinetics = TempFldr + 'Exch_Type' + str(iExch-1) + '.dat' 
            csvkinetics  = open(ExchKinetics, 'w')

            print('      [Write_Kinetics]: Writing Exchange Nb. '+ str(iExch-1) + ': ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name  + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name  )
            InelFile     = InputData.Kin.ReadFldr                               + '/'  + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name  + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name + '.csv'

            for iLevel in range(Syst.Molecule[0].NBins):
                for jLevel in range(Syst.Molecule[Syst.ExchtoMol[iExch-2]].NBins):

                    if ((Syst.T[iT-1].ProcExch[iExch-2].Rates[iLevel,jLevel] > 0.0) and (Syst.Molecule[0].LevelEEh[iLevel] > Syst.Molecule[Syst.ExchtoMol[iExch-2]].LevelEEh[jLevel]) ):
                        ProcName = Syst.Molecule[0].Name + '(' + str(iLevel+1) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '(' + str(jLevel+1) + ')+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,6\n' % Syst.T[iT-1].ProcExch[iExch-2].Rates[iLevel,jLevel]
                        csvkinetics.write(Line)

            csvkinetics.close()


def Write_Kinetics_FromOverall(Syst, Temp, InputData):
    
    mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' ) 
    mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong )
    mkdirs(    InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + '/T' + str(Temp.TranVec[iT-1]) + 'K/' )    
    TempFldr = InputData.Kin.WriteFldr + '/kinetics/' + Syst.NameLong + '/T' + str(Temp.TranVec[iT-1]) + 'K/'

    for iT in Temp.iTVec:
        print('\nTemperature Nb ', iT, '; T = ', Temp.TranVec[iT-1], 'K')

        if (InputData.Kin.WriteInel_Flg == True):
            InelKinetics = TempFldr + '/Inel.dat' 
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
                ExchKinetics = TempFldr + '/Exch_Type' + str(iExch-1) + '.dat' 
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
            if (InputData.Kin.CorrFactor != 1.0):
                print('  Writing Corrected Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
                DissKinetics = TempFldr + '/Diss_Corrected.dat' 
            else:
                print('  Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
                DissKinetics = TempFldr + '/Diss.dat' 
            csvkinetics  = open(DissKinetics, 'w')

            DissFile = InputData.Kin.ReadFldr   + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '.csv'
            with open(DissFile) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                next(readCSV)
                for row in readCSV:
                    #print(row[:])

                    if (float(row[iT]) > 0.0):
                        ProcName = Syst.Molecule[0].Name + '(' + str(row[0]) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
                        TempRate = float(row[iT])
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % (TempRate * InputData.Kin.CorrFactor)
                        csvkinetics.write(Line)
                
                csvfile.close()
                csvkinetics.close()


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
            Line = Line + '\n'
            csvQSS.write(Line)
        TempMat = np.transpose( np.expand_dims( np.concatenate( [np.array([Temp.TranVec[iT-1]], float), Syst.T[iT-1].QSS.Rate] ), axis=1 ) )
        np.savetxt(csvQSS, TempMat, delimiter=',')
    csvQSS.close()


def Write_PrefJumps(Syst, Temp, InputData, iT):

    PathToFile = InputData.FinalFldr + '/PrefJumps_Inel.csv'
    print('    [Write_PrefJumps]: Writing Jumps in File: ' + PathToFile ) 
    with open(PathToFile, 'w') as csvJumps:
        Line    = '# Level1, Level2, Level3, Level4, Level5' 
        csvJumps.write(Line)
        TempMat = Syst.T[iT-1].Proc[1].PrefJumps
        np.savetxt(csvJumps, TempMat, delimiter=',')
    csvJumps.close()

    for iProc in range(2, Syst.NProcTypes):
        PathToFile = InputData.FinalFldr + '/PrefJumps_Exch_Type' + str(iProc-1) + '.csv'
        print('    [Write_PrefJumps]: Writing Jumps in File: ' + PathToFile ) 
        with open(PathToFile, 'w') as csvJumps:
            Line    = '# Level1, Level2, Level3, Level4, Level5' 
            csvJumps.write(Line)
            TempMat = Syst.T[iT-1].ProcExch[iProc-2].PrefJumps
            np.savetxt(csvJumps, TempMat, delimiter=',')
        csvJumps.close()

