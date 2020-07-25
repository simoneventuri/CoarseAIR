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
%close all
clc

global Input Syst Temp Param Kin Rates OtherSyst OtherRates



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SPECIFYING INPUT 

%% System Inputs
Input.WORKSPACE_PATH            = '/home/venturi/WORKSPACE/'

Input.Paths.ToQCTFldr           = strcat(Input.WORKSPACE_PATH, '/CoarseAIR/O3_ALL/Test/');
Input.Paths.ToKinMainFldr       = strcat(Input.WORKSPACE_PATH, '/Air_Database/Run_0D');
Input.Paths.ToHDF5Fldr          = strcat(Input.WORKSPACE_PATH, '/Air_Database/HDF5_Database/');
Input.TranVec                   = [10000];%[1500, 2500, 5000, 6000, 8000, 10000, 12000, 14000, 15000, 20000];
Input.SystNameLong              = 'O3_UMN';
Input.iPES                      = 0;
Input.Suffix                    = ''
Input.RunSuffix                 = '';

Input.Kin.MolResolutionIn       = [{'StS'}];
Input.Kin.EqNStatesIn           = [   6115];
Input.Kin.MinStateIn            = [      1];
Input.Kin.MaxStateIn            = [   6115];
Input.Kin.PathToMappingIn       = [   {''}];
Input.Kin.PathToWriteMappingIn  = [   {''}];
Input.Kin.NGroupsIn             = [      0];
Input.Kin.MolResolutionOut      = [{'CGM'}];
Input.Kin.PathToMappingOut      = [   {''}];
Input.Kin.CGM_Strategy          = [{'CBM'}];
Input.Kin.ParamsGroupsOut       = [     2];
Input.Kin.NGroupsOut            = [     20]; %45
% Input.Kin.CGM_Strategy          = [{'RVE'}];
% Input.Kin.ParamsGroupsOut       = [     30];
% Input.Kin.NGroupsOut            = [     45]; %45
Input.Kin.PathToWriteMappingOut = [{''}];%[{'/home/venturi/WORKSPACE/Air_Database/Run_0D/database/grouping/'}];

Input.Kin.Proc.DissFlg          = 0;
Input.Kin.NBinsSuffix           = 0;
Input.Kin.DissCorrFactor        = 16.0/3.0;
Input.Kin.Proc.DissInelFlg      = 0;
Input.Kin.Proc.InelFlg          = 1;
Input.Kin.Proc.ExchFlg1         = 1;
Input.Kin.Proc.ExchFlg2         = 0;

Input.Kin.ReadRatesProc         = [2, 2, 2]
Input.Kin.RateSource            = 'HDF5'; % CoarseAIR / CG-QCT / HDF5 / PLATO
Input.Kin.ReadOtherSyst         = []
Input.Kin.OtherSystInHDF5       = []

Input.FigureFormat              = 'PrePrint';
Input.ReLoad                    = 1;


%% Inputs for Plotting
Input.iFig               = 101;
Input.SaveFigsFlgInt     = 0;
Input.Paths.SaveFigsFldr = strcat(Input.WORKSPACE_PATH, '/Air_Paper/Figures/');


%% Inputs for Saving Data
Input.Paths.SaveDataFldr = strcat(Input.WORKSPACE_PATH, '/Air_Paper/Data/');


%% Tasks Inputs
Input.Tasks.All = false

%% CoarseAIR
% Plotting Diatomic Potential
Input.Tasks.Plot_DiatPot.Flg                           = false;
Input.Tasks.Plot_DiatPot.MoleculesOI                   = [1];
Input.Tasks.Plot_DiatPot.Extremes                      = [1.5, 8.0; 1.5, 6.0];
Input.Tasks.Plot_DiatPot.jqnVec                        = [0, 100, 200];
% Plotting Overall Rate Coefficients (Dissociation and Exchange)
Input.Tasks.Plot_OverallRates.Flg                      = false;
% Plotting Pair Contributions to Dissociation Rate Coefficients
Input.Tasks.Plot_DifferentDissRates.Flg                = false;
% Writing Rates for Paraview
Input.Tasks.Write_RatesParaview.Flg                    = false;
Input.Tasks.Write_RatesParaview.MinRate                = [1e-12, 1e-12, 1.e-12]
Input.Tasks.Write_RatesParaview.Proc                   = [false, true, true]
% Input.Tasks.Write_RatesParaview.vqns                   = [0, 10,  0,20,40, 30,60, 30, 20, 20,  7, 10, 30, 25, 45, 5, 10, 10]
% Input.Tasks.Write_RatesParaview.jqns                   = [0,150,240,30,60,120,10,180,150,170,120, 50,  0, 90,110,90,210,250]
Input.Tasks.Write_RatesParaview.vqns                   = [0,  0,  0,  0,  0,  0,  8,  8,  8,  8,  8,  8, 20,20, 20, 20, 20, 30,30, 40, 40, 60,60]
Input.Tasks.Write_RatesParaview.jqns                   = [0,120,180,210,240,280, 30,100,150,180,210,240,  0,90,130,160,190, 30,90,130,160, 15,60]
Input.Tasks.Write_RatesParaview.IncludeExch            = false
% Writing Rates for Clustering
Input.Tasks.Write_RatesForClustering.Flg               = false;
Input.Tasks.Write_RatesForClustering.MinRate           = 1.e-16;
Input.Tasks.Write_RatesForClustering.WriteFldr         = strcat('/home/venturi/WORKSPACE/SpectralCluster/data/');
Input.Tasks.Write_RatesForClustering.MinState          = 1;
Input.Tasks.Write_RatesForClustering.MaxState          = 100000;
Input.Tasks.Write_RatesForClustering.IncludeExch       = true;
% Writing Rates as A Network
Input.Tasks.Write_RatesAsNetwork.Flg                   = false;
Input.Tasks.Write_RatesAsNetwork.MinRate               = 5.e-13;
Input.Tasks.Write_RatesAsNetwork.WriteFldr             = strcat('/home/venturi/WORKSPACE/SpectralCluster/data/');
Input.Tasks.Write_RatesAsNetwork.MinState              = 1;
Input.Tasks.Write_RatesAsNetwork.MaxState              = 100000;
Input.Tasks.Write_RatesAsNetwork.IncludeExch           = true;
% Compute Grouped Rate Coefficients
Input.Tasks.Compute_GroupedRates.Flg                   = false;
% Plotting Reconstructed Rate Coefficients
Input.Tasks.Plot_ReconstructedRates.Flg                = false;

%% KONIG and PLATO
% Plotting Mole Fractions
Input.Tasks.Plot_MoleFracs.Flg                         = false;
Input.Tasks.Plot_MoleFracs.CompStart                   = 1;
Input.Tasks.Plot_MoleFracs.CompEnd                     = 2;
Input.Tasks.Plot_MoleFracs.Normalize                   = 0;
% Plotting Global Rates
Input.Tasks.Plot_GlobalRates.Flg                       = false;
Input.Tasks.Plot_GlobalRates.MoleculesOI               = [1];
% Plotting Mole Fractions and Global Rates
Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg         = false;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.CompStart   = 2;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.CompEnd     = 2;
Input.Tasks.Plot_MoleFracs_and_GlobalRates.MoleculesOI = [1];
% Plotting Vib. Distribution Function
Input.Tasks.Plot_VDF.Flg                               = false;
Input.Tasks.Plot_VDF.MoleculesOI                       = [1];
Input.Tasks.Plot_VDF.tSteps                            = [1.e-14, 1e-12, 1e-10, 1e-8, 1e-6]%[8.94e-7]%[7.e-6, 30e-6, 100e-6, 5.e-3];
% Plotting RVS Populations
Input.Tasks.Plot_Populations.Flg                       = true;
Input.Tasks.Plot_Populations.MoleculesOI               = [1];
Input.Tasks.Plot_Populations.tSteps                    = [1.e-14, 1.e-13, 1e-12, 1.e-11, 1e-10, 1.e-9, 1e-8, 1e-7, 1e-6, 1e-5]%[8.94e-7]%[7.e-6, 30e-6, 100e-6, 5.e-3];
Input.Tasks.Plot_Populations.GroupColors               = 0;
Input.Tasks.Plot_Populations.ColorIdx                  = 1;
% Plotting RVS Populations Vqn Specific
Input.Tasks.Plot_PopulationsVqnSpecific.Flg            = false;
Input.Tasks.Plot_PopulationsVqnSpecific.MoleculesOI    = [1];
Input.Tasks.Plot_PopulationsVqnSpecific.tSteps         = [1.e-14, 1.e-13, 1e-12, 1.e-11, 1e-10, 1.e-9, 1e-8, 1e-7, 1e-6, 1e-5]%[8.94e-7]%[7.e-6, 30e-6, 100e-6, 5.e-3];
Input.Tasks.Plot_PopulationsVqnSpecific.GroupColors    = 2;
% Plotting Energies
Input.Tasks.Plot_Energies.Flg                          = false;
Input.Tasks.Plot_Energies.MoleculesOI                  = [1];
Input.Tasks.Plot_Energies.LTFlag                       = false;
% Plotting Energy Depletions
Input.Tasks.Plot_EnergyDepletions.Flg                  = false;
Input.Tasks.Plot_EnergyDepletions.MoleculesOI          = [1];
Input.Tasks.Plot_EnergyDepletions.RemovalProc          = [1];
Input.Tasks.Plot_EnergyDepletions.Proj                 = [1,1];
Input.Tasks.Plot_EnergyDepletions.Targ                 = [2];

% if (Input.Tasks.All)
%     Input.Kin.Proc.DissFlg          = 13;
%     Input.Kin.NBinsSuffix           = 0;
%     Input.Kin.DissCorrFactor        = 16.0/3.0;
%     Input.Kin.Proc.DissInelFlg      = 0;
%     Input.Kin.Proc.InelFlg          = 1;
%     Input.Kin.Proc.ExchFlg1         = 1;
%     Input.Kin.Proc.ExchFlg2         = 0;
% 
%     Input.Kin.ReadRatesProc         = [2, 2, 2, 0]
% 
%     Input.Tasks.Plot_Populations.Flg                       = true;
%     Input.Tasks.Plot_Energies.Flg                          = false;
% 
%     Input.Tasks.Plot_MoleFracs.Flg                         = true;
%     Input.Tasks.Plot_GlobalRates.Flg                       = true;
%     Input.Tasks.Plot_EnergyDepletions.Flg                  = true;
% else
%     Input.Kin.Proc.DissFlg          = 0;
%     Input.Kin.NBinsSuffix           = 0;
%     Input.Kin.DissCorrFactor        = 16.0/3.0;
%     Input.Kin.Proc.DissInelFlg      = 0;
%     Input.Kin.Proc.InelFlg          = 1;
%     Input.Kin.Proc.ExchFlg1         = 1;
%     Input.Kin.Proc.ExchFlg2         = 0;
%     
%     Input.Kin.ReadRatesProc         = [0, 2, 2, 0]
% 
%     Input.Tasks.Plot_Populations.Flg                       = false;
%     Input.Tasks.Plot_Energies.Flg                          = false;
% 
%     Input.Tasks.Plot_MoleFracs.Flg                         = false;
%     Input.Tasks.Plot_GlobalRates.Flg                       = false;
%     Input.Tasks.Plot_EnergyDepletions.Flg                  = false;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initializing
Syst.NameLong = Input.SystNameLong;
Syst          = Initialize_ChemicalSyst(Syst)
Initialize_Parameters()
Initialize_Input()


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
pause

iFigStart = Input.iFig;
%% Looping On Translational Temperatures
for iT = 1:length(Temp.TranVec)
    Temp.iT       = iT;
    Temp.TNow     = Temp.TranVec(iT);
    Temp.TNowChar = num2str(Temp.TranVec(iT));
    Input.iFig    = iFigStart;
    Input.Paths.ToKinRunFldr = strcat(Input.Paths.ToKinMainFldr, '/output_', Syst.NameLong, Input.RunSuffix, '_T', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg);
    
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
            Input.Tasks.Write_RatesForClustering.Flg       || ...
            Input.Tasks.Write_RatesAsNetwork.Flg           || ...
            Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_Populations.Flg               || ...
            Input.Tasks.Plot_PopulationsVqnSpecific.Flg    || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_Energies.Flg                  || ...
            Input.Tasks.Plot_EnergyDepletions.Flg          || ...
            (sum(Input.Kin.ReadRatesProc)>0) ) 
        
            %% Reading Rates
            Read_Rates()
            
        end
        
        if (Input.Tasks.Plot_MoleFracs.Flg                 || ...
            Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_VDF.Flg                       || ...
            Input.Tasks.Plot_Populations.Flg               || ...
            Input.Tasks.Plot_PopulationsVqnSpecific.Flg    || ...
            Input.Tasks.Plot_Energies.Flg                  || ...
            Input.Tasks.Plot_EnergyDepletions.Flg)
        
            %% Reading Thermodynamics Variables Outputted by KONIG
            Read_KONIGBox() 
            
        end
        
        if (Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_VDF.Flg                       || ...
            Input.Tasks.Plot_Populations.Flg               || ...
            Input.Tasks.Plot_PopulationsVqnSpecific.Flg    || ...
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
        
        if ( (sum(Input.Kin.ReadRatesProc)>0) )
        
            %% Computing Thermal Rates
            Compute_Rates_Thermal()   
        
        end
        
        if ((Input.Tasks.Plot_Populations.Flg            && Input.Kin.Proc.DissFlg > 0) || ...
            (Input.Tasks.Plot_PopulationsVqnSpecific.Flg && Input.Kin.Proc.DissFlg > 0) || ...
            Input.Tasks.Plot_GlobalRates.Flg               || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg || ...
            Input.Tasks.Plot_EnergyDepletions.Flg)

            %% Computing Thermal Rates
            Compute_Rates_Global()   
        
        end
        
        if ((Input.Tasks.Plot_Populations.Flg            && Input.Kin.Proc.DissFlg > 0) || ...
            (Input.Tasks.Plot_PopulationsVqnSpecific.Flg && Input.Kin.Proc.DissFlg > 0) || ...
            Input.Tasks.Plot_GlobalRates.Flg                                 || ...
            Input.Tasks.Plot_MoleFracs_and_GlobalRates.Flg                   || ...
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
    
    
    %%%% Writing Quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  
    
    %% Writing Rate Coefficients for Paraview
    if (Input.Tasks.Write_RatesParaview.Flg)
        Write_RatesForParaview(Input.Tasks.Write_RatesParaview)
    end
    
    %% Writing Rate Coefficients for Paraview
    if (Input.Tasks.Write_RatesAsNetwork.Flg)
        Write_RatesAsNetwork(Input.Tasks.Write_RatesAsNetwork)
    end
        
    %% Writing Rate Coefficients for Clustering
    if (Input.Tasks.Write_RatesForClustering.Flg)
        Write_RatesForClustering(Input.Tasks.Write_RatesForClustering)
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
    
    %% Plotting RVS Populations
    if (Input.Tasks.Plot_PopulationsVqnSpecific.Flg)
       Plot_PopulationsVqnSpecific(Input.Tasks.Plot_PopulationsVqnSpecific) 
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
    clear Rates Kin
end