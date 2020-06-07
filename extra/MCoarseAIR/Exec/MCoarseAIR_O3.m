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

global Input Syst Temp Param Kin Rates OtherSyst OtherRates



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SPECIFYING INPUT 

%% System Inputs
Input.WORKSPACE_PATH        = '/home/venturi/WORKSPACE/'

Input.Paths.ToQCTFldr       = strcat(Input.WORKSPACE_PATH, '/CoarseAIR/O3_ALL/Test/');
Input.Paths.ToKinMainFldr   = strcat(Input.WORKSPACE_PATH, '/O3Diss_Database/Run_0D');
Input.Paths.ToHDF5Fldr      = strcat(Input.WORKSPACE_PATH, '/O3Diss_Database/HDF5_Database/');
Input.TranVec               = [5000, 10000, 20000]%[5000, 6000, 8000, 10000, 12000, 14000, 15000, 20000];
Input.SystNameLong          = 'O3_UMN';
Input.iPES                  = 0;
Input.Suffix                = ''
Input.RunSuffix             = '';

Input.Kin.MolResolutionIn   = ['StS'];
Input.Kin.EqNStatesIn       = [ 6115];
Input.Kin.MinStateIn        = [    1];
Input.Kin.MaxStateIn        = [ 6115];
Input.Kin.PathToMappingIn   = [   ''];
Input.Kin.NGroupsIn         = [    0];
Input.Kin.MolResolutionOut  = ['CGM'];
Input.Kin.PathToMappingOut  = [   ''];
Input.Kin.CGM_Strategy      = ['CBM'];
Input.Kin.ParamsGroupsOut   = [  1.0];
Input.Kin.NGroupsOut        = [   45];

Input.Kin.Proc.DissFlg      = 2;
Input.Kin.NBinsSuffix       = 0;
Input.Kin.DissCorrFactor    = 16.0/3.0;
Input.Kin.Proc.DissInelFlg  = 0;
Input.Kin.Proc.InelFlg      = 1;
Input.Kin.Proc.ExchFlg1     = 0;
Input.Kin.Proc.ExchFlg2     = 0;

Input.Kin.ReadRatesProc     = [false, false, false]
Input.Kin.RateSource        = 'HDF5'; % CoarseAIR / CG-QCT / HDF5 / PLATO
Input.Kin.ReadOtherSyst     = []
Input.Kin.OtherSystInHDF5   = []

Input.FigureFormat          = 'PrePrint';
Input.ReLoad                = 1;


%% Inputs for Plotting
Input.iFig               = 101;
Input.SaveFigsFlgInt     = 0;
Input.Paths.SaveFigsFldr = strcat(Input.WORKSPACE_PATH, '/Mars_Paper/Figures/');


%% Inputs for Saving Data
Input.Paths.SaveDataFldr = strcat(Input.WORKSPACE_PATH, '/Mars_Paper/Data/');


%% Tasks Inputs

%% CoarseAIR
% Plotting Diatomic Potential
Input.Tasks.Plot_DiatPot.Flg                           = true;
Input.Tasks.Plot_DiatPot.MoleculesOI                   = [1];
Input.Tasks.Plot_DiatPot.Extremes                      = [1.5, 8.0; 1.5, 6.0];
Input.Tasks.Plot_DiatPot.jqnVec                        = [0, 100, 200];
% Plotting Overall Rate Coefficients (Dissociation and Exchange)
Input.Tasks.Plot_OverallRates.Flg                      = false;
% Plotting Pair Contributions to Dissociation Rate Coefficients
Input.Tasks.Plot_DifferentDissRates.Flg                = false;
% Writing Rates for Paraview
Input.Tasks.Write_RatesParaview.Flg                    = false;
Input.Tasks.Write_RatesParaview.MinRate                = [1e-12, 1e-12, 1e-12]
% Compute Grouped Rate Coefficients
Input.Tasks.Compute_GroupedRates.Flg                   = false;
% Plotting Reconstructed Rate Coefficients
Input.Tasks.Plot_ReconstructedRates.Flg                = false;

%% KONIG and PLATO
% Plotting Mole Fractions
Input.Tasks.Plot_MoleFracs.Flg                         = false;
Input.Tasks.Plot_MoleFracs.CompStart                   = 1;
Input.Tasks.Plot_MoleFracs.CompEnd                     = 2;
% Plotting Global Rates
Input.Tasks.Plot_GlobalRates.Flg                       = false;
Input.Tasks.Plot_GlobalRates.MoleculesOI               = [1,2];
% Plotting Mole Fractions and Global Rates
Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg         = false;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.CompStart   = 2;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.CompEnd     = 2;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.MoleculesOI = [1,2];
% Plotting Vib. Distribution Function
Input.Tasks.Plot_VDF.Flg                               = false;
Input.Tasks.Plot_VDF.MoleculesOI                       = [1];
Input.Tasks.Plot_VDF.tSteps                            = [8.0e-7]%[1.23e-6]%[8.94e-7]%[7.e-6, 30e-6, 100e-6, 5.e-3];
% Plotting RVS Populations
Input.Tasks.Plot_Populations.Flg                       = false;
Input.Tasks.Plot_Populations.MoleculesOI               = [1];
Input.Tasks.Plot_Populations.tSteps                    = [1.23e-6]%[8.94e-7]%[7.e-6, 30e-6, 100e-6, 5.e-3];
Input.Tasks.Plot_Populations.GroupColors               = 2;
% Plotting Energies
Input.Tasks.Plot_Energies.Flg                          = false;
Input.Tasks.Plot_Energies.MoleculesOI                  = [1];
Input.Tasks.Plot_Energies.LTFlag                       = true;
% Plotting Energy Depletions
Input.Tasks.Plot_EnergyDepletions.Flg                  = false;
Input.Tasks.Plot_EnergyDepletions.MoleculesOI          = [1];
Input.Tasks.Plot_EnergyDepletions.RemovalProc          = [1];
Input.Tasks.Plot_EnergyDepletions.ProjTarg             = [1,2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing
Syst.NameLong = Input.SystNameLong;
Syst          = Initialize_ChemicalSyst(Syst)
Initialize_Input()
Initialize_Parameters()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reading Quantities
if Input.ReLoad > 0 

    %% Reading Levels Info
    Syst = Read_LevelInfo(Syst)
    for iSyst = 1:length(Input.Kin.ReadOtherSyst)
        if (Input.Kin.ReadOtherSyst(iSyst))
            OtherSyst(iSyst).Syst = Read_LevelInfo(OtherSyst(iSyst).Syst);
        end
    end

    %% Grouping the Levels in Output
    Group_Out()

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
        %close all
        
        
        %%%% Reading Quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        
        %% Reading Group Energies and Part Funcs
        Syst = Read_EeV_and_Q_CG(Syst) 
        for iSyst = 1:length(Input.Kin.ReadOtherSyst)
            if (Input.Kin.ReadOtherSyst(iSyst))
                OtherSyst(iSyst).Syst = Read_EeV_and_Q_CG(OtherSyst(iSyst).Syst);
            end
        end
        
        %% Compute Equilibrium Constants
        Compute_EqConsts()
        
        if (Input.Tasks.Plot_OverallRates.Flg              || ...
            Input.Tasks.Plot_DifferentDissRates.Flg        || ...
            Input.Tasks.Write_RatesParaview.Flg            || ...
            Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_Energies.Flg                  || ...
            Input.Tasks.Plot_EnergyDepletions.Flg          || ...
            any(Input.Kin.ReadRatesProc))
        
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
        
        if (Input.Tasks.Compute_GroupedRates.Flg)
           
           %% Grouping the Rate Coefficients
           Compute_GroupedRates()
            
        end
        
        if (any(Input.Kin.ReadRatesProc(1,:)))
        
            %% Computing Thermal Rates
            Compute_Rates_Thermal()   
        
        end
        
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
    
    %% Plotting Reconstructed Rate Coefficients
    if (Input.Tasks.Plot_ReconstructedRates.Flg)
        Plot_ReconstructedRates()
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
    
    
    %pause
    %clear Rates Kin
end