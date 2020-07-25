%close all
clear all
clc

global Input Syst Params

Input.SystNameLong = 'N4_NASA';
Input.FigureFormat = 'PrePrint';

Syst.NameLong = Input.SystNameLong;
Syst          = Initialize_ChemicalSyst(Syst)
Initialize_Parameters()


R       = linspace(1.5,10.0, 3000);
[V_UMN, dV]  = N2_UMN_PIPNN_ForN4(R); 
hold on
[V_NASA, dV] = N2_NASA(R); 


figure
plot(R,V_UMN)
hold on
plot(R,V_NASA)
