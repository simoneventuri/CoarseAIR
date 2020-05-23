% -- MATLAB --
%%==============================================================================================================
% 
% Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
% 
% Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
%
% Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
% 
% This program is free software; you can redistribute it and/or modify it under the terms of the 
% Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without e=ven the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% See the GNU Lesser General Public License for more details. 
% 
% You should have received a copy of the GNU Lesser General Public License along with this library; 
% if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
% 
%---------------------------------------------------------------------------------------------------------------
%%==============================================================================================================

clear all
close all
clc

global Input Syst Temp Param Kin Rates



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SPECIFYING INPUT 

%% System Inputs
Input.Paths.ToQCTFldr       = '/home/venturi/WORKSPACE/CoarseAIR/O2C_ALL/Test/';
Input.Paths.ToKinMainFldr   = '/home/venturi/WORKSPACE/Mars_Database/Run_0D/';
Input.Paths.ToHDF5Fldr      = '/home/venturi/WORKSPACE/Mars_Database/HDF5_Database/';
Input.TranVec               = 20000 %[2500 5000 7500 10000 12500 15000 20000];
Input.SystNameLong          = 'O2C_NASA';
Input.iPES                  = 0;
Input.Suffix                = ''
Input.Kin.MolResolutionIn   = ['StS'; 'StS'];
Input.Kin.MinStateIn        = [    1,     1];
Input.Kin.MaxStateIn        = [ 6078; 13521];
Input.Kin.NGroupsIn         = [    0,     0];
Input.Kin.MolResolutionOut  = ['StS'; 'StS'];
Input.Kin.PathToMappingOut  = [''];
Input.Kin.Proc.DissFlg      = 1;
Input.Kin.NBinsSuffix       = 0;
Input.Kin.DissCorrFactor    = 1.0;
Input.Kin.Proc.DissInelFlg  = 0;
Input.Kin.Proc.InelFlg      = 1;
Input.Kin.Proc.ExchFlg1     = 0;
Input.Kin.Proc.ExchFlg2     = 0;
Input.Kin.RateSource        = 'HDF5'; % CoarseAIR / CG-QCT / HDF5 / PLATO
Input.FigureFormat          = 'PrePrint';
Input.ReLoad                = 1;

Input.RunSuffix = '';


%% Inputs for Plotting
Input.iFig               = 101;
Input.SaveFigsFlgInt     = 0;
Input.Paths.SaveFigsFldr = '/home/venturi/WORKSPACE/Mars_Paper/Figures/temp/';


%% Inputs for Saving Data
Input.Paths.SaveDataFldr = '/home/venturi/WORKSPACE/Mars_Paper/Data/temp/';


%% Tasks Inputs

%% CoarseAIR
% Plotting Diatomic Potential
Input.Tasks.Plot_DiatPot.Flg                           = false;
Input.Tasks.Plot_DiatPot.Extremes                      = [1.5, 10.0; 1.5, 10.0];
Input.Tasks.Plot_DiatPot.jqnVec                        = [44];
% Plotting Overall Rate Coefficients (Dissociation and Exchange)
Input.Tasks.Plot_OverallRates.Flg                      = false;
% Plotting Pair Contributions to Dissociation Rate Coefficients
Input.Tasks.Plot_DifferentDissRates.Flg                = false;
% Writing Rates for Paraview
Input.Tasks.Write_RatesParaview.Flg                    = false;


%% KONIG and PLATO
% Plotting Mole Fractions
Input.Tasks.Plot_MoleFracs.Flg                         = false;
Input.Tasks.Plot_MoleFracs.CompStart                   = 1;
Input.Tasks.Plot_MoleFracs.CompEnd                     = 4;
% Plotting Global Rates
Input.Tasks.Plot_GlobalRates.Flg                       = false;
% Plotting Mole Fractions and Global Rates
Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg         = true;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.CompStart   = 1;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.CompEnd     = 4;
% Plotting Vib. Distribution Function
Input.Tasks.Plot_VDF.Flg                               = false;
Input.Tasks.Plot_VDF.MoleculesOI                       = [1];
Input.Tasks.Plot_VDF.tSteps                            = [8.0e-7]%[1.23e-6]%[8.94e-7]%[7.e-6, 30e-6, 100e-6, 5.e-3];
% Plotting RVS Populations
Input.Tasks.Plot_Populations.Flg                       = true;
Input.Tasks.Plot_Populations.MoleculesOI               = [1];
Input.Tasks.Plot_Populations.tSteps                    = [1.e-10, 8.e-7, 1.e-6];
Input.Tasks.Plot_Populations.GroupColors               = 2;
% Plotting Energies
Input.Tasks.Plot_Energies.Flg                          = false;
Input.Tasks.Plot_Energies.MoleculesOI                  = [1];
Input.Tasks.Plot_Energies.LTFlag                       = true;
% Plotting Energy Depletions
Input.Tasks.Plot_EnergyDepletions.Flg                  = true;
Input.Tasks.Plot_EnergyDepletions.MoleculesOI          = [1];
Input.Tasks.Plot_EnergyDepletions.RemovalProc          = [1];
Input.Tasks.Plot_EnergyDepletions.ProjTarg             = [2,3];

Input.Tasks.ReadRates          = true;
Input.Tasks.ReadInelasticRates = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing
Initialize_ChemicalSyst()
Initialize_Input()
Initialize_Parameters()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing
Initialize_ChemicalSyst()
Initialize_Input()
Initialize_Parameters()



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reading Quantities
if Input.ReLoad > 0 

    %% Reading Levels Info
    Read_LevelInfo()

end


iFigStart = Input.iFig;
%% Looping On Translational Temperatures
for iT = 1:length(Temp.TranVec)
    Temp.iT       = iT;
    Temp.TNow     = Temp.TranVec(iT);
    Temp.TNowChar = num2str(Temp.TranVec(iT));
    Input.iFig    = iFigStart;

    Input.Paths.ToKinRunFldr = strcat(Input.Paths.ToKinMainFldr, '/output_', Syst.NameLong, '_T', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg, Input.RunSuffix);
    
    
    if Input.ReLoad > 0 
        
        
        %%%% Reading Quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        
        %% Reading Group Energies and Part Funcs
        Read_EeV_and_Q_CG() 
        
        if (Input.Tasks.Plot_OverallRates.Flg              || ...
            Input.Tasks.Plot_DifferentDissRates.Flg        || ...
            Input.Tasks.Write_RatesParaview.Flg            || ...
            Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_Energies.Flg                  || ...
            Input.Tasks.Plot_EnergyDepletions.Flg          || ...
            Input.Tasks.ReadRates)
        
            %% Reading Rates
            Read_Rates()
            
        end
        
        if (Input.Tasks.Plot_MoleFracs.Flg                 || ...
            Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_VDF.Flg                       || ...
            Input.Tasks.Plot_Populations.Flg               || ...
            Input.Tasks.Plot_Energies.Flg                  || ...
            Input.Tasks.Plot_EnergyDepletions.Flg)
        
            %% Reading Thermodynamics Variables Outputted by KONIG
            Read_KONIGBox() 
            
        end
        
        if (Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_VDF.Flg                       || ...
            Input.Tasks.Plot_Populations.Flg               || ...
            Input.Tasks.Plot_Energies.Flg                  || ...
            Input.Tasks.Plot_EnergyDepletions.Flg)
        
            %% Reading Level/Group Population Outputted by KONIG
            Read_Pops()    
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%% Computing Quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        
        %% Computing Thermal Rates
        %Compute_Rates_Thermal()   
        
        if (Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg)

            %% Computing Thermal Rates
            Compute_Rates_Global()   
        
        end
        
        if (Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_EnergyDepletions.Flg)

            %% Computing Rate Values and Initial-Final Times for QSS 
            Compute_QSS()
            
        end
        
        if (Input.Tasks.Plot_Energies.Flg                  || ...
            Input.Tasks.Plot_EnergyDepletions.Flg)
        
            %% Computing Energies
            Compute_Energies(Input.Tasks.Plot_EnergyDepletions)
        
        end
        
        if (Input.Tasks.Plot_EnergyDepletions.Flg)
            
            %% Computing Energy Depletions
            Compute_EnergyDepletions(Input.Tasks.Plot_EnergyDepletions)
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    
    %%%% Plotting Quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   
    
    %% Plotting Diatomic Potential
    if (Input.Tasks.Plot_DiatPot.Flg)
        Plot_DiatPot(Input.Tasks.Plot_DiatPot)
    end
    
    %% Plotting Overall Rate Coefficients (Dissociation and Exchange)
    if (Input.Tasks.Plot_OverallRates.Flg)
        Plot_OverallRates()    
    end
    
    %% Plotting Pair Contributions to Dissociation Rate Coefficients
    if (Input.Tasks.Plot_DifferentDissRates.Flg)
        Plot_DifferentDissRates()
    end
    
    %% Writing Rate Coefficients for Paraview
    if (Input.Tasks.Write_RatesParaview.Flg)
        Write_RatesForParaview(Input.Tasks.Write_RatesParaview)
    end
    
    
    %% Plotting Mole Fractions
    if (Input.Tasks.Plot_MoleFracs.Flg)
        Plot_MoleFracs(Input.Tasks.Plot_MoleFracs)
    end
    
    %% Plotting Global Rates (Dissociation and Exchange)
    if (Input.Tasks.Plot_GlobalRates.Flg)
        Plot_GlobalRates(Input.Tasks.Plot_GlobalRates)    
    end
    
    %% Plotting Global Rates (Dissociation and Exchange) on top of Mole Fractions
    if (Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg)
        Plot_MoleFracs_and_GlobalRates(Input.Tasks.Plot_MoleFracs_and_GlobalRates)
    end
    
    %% Plotting Vib. Distr. Function
    if (Input.Tasks.Plot_VDF.Flg)
       Plot_VDF(Input.Tasks.Plot_VDF) 
    end
        
    %% Plotting RVS Populations
    if (Input.Tasks.Plot_Populations.Flg)
       Plot_Populations(Input.Tasks.Plot_Populations) 
    end
    
    %% Plotting Energies
    if (Input.Tasks.Plot_Energies.Flg)
        Plot_Energies(Input.Tasks.Plot_Energies)
    end
    
    %% Plotting Energy Depletions
    if (Input.Tasks.Plot_EnergyDepletions.Flg)
        Plot_EnergyDepletions(Input.Tasks.Plot_EnergyDepletions)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pause
    clear Rates Kin
end