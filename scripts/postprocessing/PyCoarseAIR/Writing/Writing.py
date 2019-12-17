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

def Write_PartFuncsAndEnergies(Syst, Temp, InputData):

    for iMol in range(Syst.NMolecules):

        for iT in Temp.iTVec:

            PathToFileOrig = InputData.WriteKinFolder + '/thermo/' + Syst.Molecule[iMol].Name + '_Format'
            PathToFile     = InputData.WriteKinFolder + '/thermo/' + Syst.Molecule[iMol].Name + '_' + str(Temp.TranVec[iT-1])
            DestTemp       = shutil.copyfile(PathToFileOrig, PathToFile)
            print('Copied File to: ', DestTemp)

            with open(PathToFile, 'a') as f:
                Line = 'NB_ENERGY_LEVELS = ' + str(Syst.Molecule[iMol].NBins) + '\n'
                f.write(Line)
                np.savetxt(f, np.transpose(np.array([Syst.Molecule[iMol].T[iT-1].Q, Syst.Molecule[iMol].LevelEeV])), fmt='%.8e    %.8e')



def Write_Kinetics_FromOverall(Syst, Temp, InputData):

    for iT in Temp.iTVec:
        print('\nTemperature Nb ', iT, '; T = ', Temp.TranVec[iT-1], 'K')

        if (InputData.WriteInelKin_Flg == True):
            InelKinetics = InputData.WriteKinFolder + '/kinetics/' + Syst.Name + 'Inel_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
            csvkinetics  = open(InelKinetics, 'w')

            print('  Writing Inelastic: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name )
            InelFile     = InputData.ReadKinFolder  + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '.csv'
            with open(InelFile) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                next(readCSV)
                for row in readCSV:

                    if (float(row[iT+1]) > 0.0):
                        ProcName = Syst.Molecule[0].Name + '(' + str(row[0]) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[0].Name + '(' + str(row[1]) + ')+' + Syst.Atom[2].Name
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,5\n' % float(row[iT+1])
                        csvkinetics.write(Line)
            
                csvkinetics.close()


        if (InputData.WriteExchKin_Flg == True):

            for iExch in range (2, Syst.NProcTypes):
                print('  Writing Exchange: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name )
                ExchKinetics = InputData.WriteKinFolder + '/kinetics/' + Syst.Name + 'Exch_Type' + str(iExch-1) + '_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
                csvkinetics  = open(ExchKinetics, 'w')

                ExchFile     = InputData.ReadKinFolder + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name + '_Exch.csv'
                with open(ExchFile) as csvfile:
                    readCSV = csv.reader(csvfile, delimiter=',')
                    next(readCSV)
                    for row in readCSV:

                        if (float(row[iT+1]) > 0.0):
                            ProcName = Syst.Molecule[0].Name + '(' + str(row[0]) + ')+' + Syst.Atom[2].Name + '=' + Syst.Molecule[Syst.ExchtoMol[iExch-2]].Name + '(' + str(row[1]) + ')+' + Syst.Atom[Syst.ExchtoAtom[iExch-2]].Name
                            Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,6\n' % float(row[iT+1])
                            csvkinetics.write(Line)
                
                    csvkinetics.close()


        if (InputData.WriteDissKin_Flg == True):
            DissKinetics = InputData.WriteKinFolder + '/kinetics/' + Syst.Name + 'Diss_' + str(Temp.TranVec[iT-1]) + 'K.dat' 
            csvkinetics  = open(DissKinetics, 'w')

            print('  Writing Dissociation: ' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name )
            DissFile = InputData.ReadKinFolder   + '/' + Syst.Molecule[0].Name + '+' + Syst.Atom[2].Name + '_' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name + '.csv'
            with open(DissFile) as csvfile:
                readCSV = csv.reader(csvfile, delimiter=',')
                next(readCSV)
                for row in readCSV:
                    #print(row[:])

                    if (float(row[iT]) > 0.0):
                        ProcName = Syst.Molecule[0].Name + '(' + str(row[0]) + ')+' + Syst.Atom[2].Name + '=' + Syst.Atom[0].Name + '+' + Syst.Atom[1].Name + '+' + Syst.Atom[2].Name
                        Line     = ProcName + ':%.4e,+0.0000E+00,+0.0000E+00,2\n' % float(row[iT])
                        csvkinetics.write(Line)
            
                csvkinetics.close()