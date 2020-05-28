close all
clear all
clc


global Param Input Syst
iFig = 1;

% Input.SystNameLong       = 'CO2_NASA'
% Input.Paths.KGlobal      = '/home/venturi/WORKSPACE/Mars_Paper/Data/CO2_NASA/KGlobal_1_1_0_0.csv'
Input.SystNameLong       = 'O2C_NASA'
Input.Paths.KGlobal      = '/home/venturi/WORKSPACE/Mars_Paper/Data/O2C_NASA/KGlobal_1_1_1_0.csv'

Input.FigureFormat       = 'PrePrint';

Input.iFig               = 101;
Input.SaveFigsFlgInt     = 0;
Input.Paths.SaveFigsFldr = strcat('/home/venturi/WORKSPACE/Mars_Paper/Figures/', Input.SystNameLong);


Initialize_ChemicalSyst()
Initialize_Parameters()


Plot_Ks()