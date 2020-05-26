close all
clear all
clc


global Param Input
iFig = 1;

Input.SystNameLong       = 'CO2_NASA'
Input.Paths.TausWOExch   = '/home/venturi/WORKSPACE/Mars_Paper/Data/CO2_NASA/CO2_NASA/Taus_CO_0_1_0_0.csv'
Input.Paths.TausWExch    = '/home/venturi/WORKSPACE/Mars_Paper/Data/CO2_NASA/CO2_NASA/Taus_CO_0_1_1_0.csv'
% Input.SystNameLong       = 'O2C_NASA'
% Input.Paths.TausWOExch   = '/home/venturi/WORKSPACE/Mars_Paper/Data/O2C_NASA/O2C_NASA/Taus_O2_0_1_0_0.csv'
% Input.Paths.TausWExch    = ''

Input.FigureFormat       = 'PrePrint';

Input.iFig               = 101;
Input.SaveFigsFlgInt     = 2;
Input.Paths.SaveFigsFldr = strcat('/home/venturi/WORKSPACE/Mars_Paper/Figures/', Input.SystNameLong);


Initialize_ChemicalSyst()
Initialize_Parameters()


Plot_Taus()