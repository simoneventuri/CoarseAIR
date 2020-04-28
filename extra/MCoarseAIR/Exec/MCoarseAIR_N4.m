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

global Input Syst Temp Param Rates Kin 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SPECIFYING INPUT 

%% System Inputs
Input.Paths.ToQCTFldr       = '/home/venturi/WORKSPACE/CoarseAIR/N4_VS_NASA/Test'
Input.Paths.ToKinMainFldr   = '/home/venturi/WORKSPACE/Mars_Database/Run_0D/'
Input.Paths.ToHDF5Fldr      = '/home/venturi/WORKSPACE/Mars_Database/HDF5_Database/'
Input.TranVec               = [20000]
Input.SystNameLong          = 'N4_NASA'
Input.Kin.MolResolutionIn   = ['VSM'] 
Input.Kin.MinStateIn        = [    1]
Input.Kin.MaxStateIn        = [   61]
Input.Kin.NGroupsIn         = [    0]
Input.Kin.Proc.DissFlg      = 0
Input.Kin.Proc.DissInelFlg  = 1
Input.Kin.Proc.InelFlg      = 1
Input.Kin.Proc.ExchFlg1     = 1
Input.Kin.Proc.ExchFlg2     = 0
Input.Kin.RateSource        = 'CoarseAIR' % CoarseAIR / HDF5 / PLATO
Input.FigureFormat          = 'PrePrint'
Input.ReLoad                = 1



%% Tasks Inputs

% Plotting Mole Fractions
Input.Tasks.Plot_MoleFracs.Flg       = true
Input.Tasks.Plot_MoleFracs.CompStart = 1
Input.Tasks.Plot_MoleFracs.CompEnd   = 1



%% Plots Inputs
Input.iFig           = 1
Input.SaveFigsFlgInt = 0



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



%% Looping On Translational Temperatures
for iT = 1:length(Temp.TranVec)
    Temp.iT       = iT;
    Temp.TNow     = Temp.TranVec(iT);
    Temp.TNowChar = num2str(Temp.TranVec(iT));
  
    Input.Paths.ToKinRunFldr = strcat(Input.Paths.ToKinMainFldr, '/output_', Syst.NameLong, '_T', Temp.TNowChar, 'K_', Input.Kin.Proc.OverallFlg);


    if Input.ReLoad > 0 
        
        %% Reading Group Energies and Part Funcs
        Read_EeV_and_Q_CG(); 
        
        %% Reading Rates
        Read_Rates();
        
        %% Reading Thermodynamics Variables Outputted by KONIG
        Read_KONIGBox() 
        
        %% Reading Level/Group Population Outputted by KONIG
        Read_Pops()    
    end
    
    
    %% Computing Thermal Rates
    Compute_Rates_Thermal()   
    
    
    %% Plotting Mole Fractions
    if (Input.Tasks.Plot_MoleFracs.Flg)
        Plot_MoleFracs(Input.Tasks.Plot_MoleFracs)
    end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% CHOOSING POSTPROCESSING FUNCTIONS
% PlotMolFracsFlg            = 1
% PlotTemperaturesFlg        = 0
% PlotKsFlg                  = 0
% PlotKsMoleFracsFlg         = 0
% PlotMovieFlg               = 0
% PlotStepsFlg               = 0
% PlotEnergiesFlg            = 0
% PlotEnergyDepletionsFlg    = 0
% PlotQSSDistrFlg            = 0
% PlotProcProbFlg            = 0
% PlotNewSahaFlg             = 0
% PlotTausFlg                = 0
% LevToBinFlg                = 0
% WriteNetworkCSVFlg         = 0
% WriteRatesToParaviewCSVFlg = 0
% ComputeMasterEquationFlg   = 0
%   NStepsMovie           = 1000;
%   TimeMinMovie(1:3)     = 1.d-10;
%   TimeMaxMovie(1)       = 1.e-6 %2.3E-07;
%   TimeMaxMovie(2)       = 1.e-6; %2.3E-07;
%   TimeMaxMovie(3)       = 1.e-6; %2.3E-07;
%   StartBin              = 1
%   FinalBin              = 6115
%   NSubplots             = 3
%   NInstantsPerSubplot   = 3
%   tInstants             = [0.0, 1.e-13, 1.e-11, 1.e-9, 5.e-9, 1.e-8, 5.e-8, 1.e-7]%[9.311e-7]%[0.0, 1.e-13, 1.e-11, 1.e-9, 5.e-9, 1.e-8, 5.e-8, 1.e-7, 5.e-7, 8.e-7, 1.e-6]%, 3.e-7]
%   MinEvPlot             = [-20.d0, -1.d2];
%   MaxEvPlot             = [1.d2,  1.d2];
%   StepsOverlappingSteps = 0
%   vqnColor              = 0
%   SubplotsFlg           = 0
%   QNSpace               = 0
%   ReadAllRatesFlg       = -1
%   SaveRatesFlg          = 1
%   XLimPlot              = [1.d-12, 1.d-1]
%   YLimPlot              = [ 0.d0,  1.d0]
%   CompOI                = 1
%   PlotPairUser          = [0, 1, 1]
%   BinnedFlg             = 5
%   MovieFlg              = 0
%   ColPartToComp         = 1
%   iFigureStart          = 2001  
%   WriteRatesFlag        = 1
%   WriteSrcTermsFlag     = 0
%   WriteDissFlag         = 0
%   WriteInelasticFlag    = 1
%   WriteAllFlag          = 0
%   
%   DissCorrectionFactor  = 1.0;
% %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ls
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% UPDATE PARAMETERS & PATHS
% [iFigure]   = UpdateParameters();
% 
% RatesPath   = strcat(PathToOutput,System(1,:),'/',MoleculesName(1,:),'/Rates/');                                                             % For Levels / Bins Rates
% SystemPath  = strcat(PathToOutput,System(1,:));   
% 
% iFigure    = 6001;
% SaveFigs   = 0
% FigDirPath = strcat('./CGQCT-',System,'-Figures/');
% if SaveFigs > 0
%   [status,msg,msgID] = mkdir(FigDirPath)
% end   
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% COMPUTING SYSTEM MAIN QUANTITIES
% ComputeSystemQuantities() 
% if ReUpload > 0
%   Rates               = zeros(TotDim(end),NBins(1),NTint);
% %   RatesSigma          = zeros(TotDim(end),NBins(1),NTint);
%   RatesMatrix         = zeros(NBins(1),max(yDim),3,NTint); 
% %   RatesSigmaMatrix    = zeros(NBins(1),max(yDim),3,NTint); 
%   DissRates           = zeros(NBins(1),NTint);
% %   DissRatesSigma      = zeros(NBins(1),NTint);
%   ProcessesRates      = zeros(NBins(1),4,NTint); 
% %   ProcessesRatesSigma = zeros(NBins(1),4,NTint); 
%   ProcessesRatesOverall = [];
% end
% 
%     
% for iT = 1:length(T0_Vec)
%     
%   if DissFlg<4 || DissFlg==6
%     OutputPath  = strcat('/home/venturi/WORKSPACE/Mars_Database/Run_0D/output_NaNbNcNd_UMN_VS_T',num2str(T0_Vec(iT)),'K_',num2str(DissFlg),'_',num2str(DissInelFlg),'_',num2str(InelFlg),'_',num2str(ExchFlg),'_',num2str(ExchFlg)); TypeCG = 0; 
%   else
%     OutputPath  = strcat('/home/venturi/WORKSPACE/Mars_Database/Run_0D/output_NaNbNcNd_UMN_VS_T',num2str(T0_Vec(iT)),'K_',num2str(DissFlg),'_',num2str(DissInelFlg),'_',num2str(InelFlg),'_',num2str(ExchFlg),'_',num2str(ExchFlg),'_',num2str(TempBins),'Bins/'); TypeCG = 0; 
%   end
%   %OutputPath  = strcat('/home/venturi/WORKSPACE/Mars_Database/Run_0D/output_O3_UMN_T20000K_3_1_1_1/'); TypeCG = 1;
%   %OutputPath  = strcat('/home/venturi/WORKSPACE/Mars_Database/Run_0D/output_O3_UMN_T20000K_4_1_1_1_45Bins/'); TypeCG = 2;
%   %OutputPath  = strcat('/home/venturi/WORKSPACE/Mars_Database/Run_0D/output_O3_UMN_T15000K_5_1_1_1_45Bins/')
% 
%   if ReUpload > 0 
% 
%     [QBins, EeV] = ReadPartFunEnergy(iT); 
%     %if PlotProcProbFlg == 1 || PlotEnergiesFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1
%       [NLevels, Levelvqn, Leveljqn, LevelEh, LevelEeV, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelg, LevelToBin, LevelQ, vToLevel, DeltaEintDiss, DeltaEintDepth, rIn, rOut, rMin, rMax, VMin, VMax, Tau, Egam,  dVIn, ddVIn, dVOut, ddVOut, dVdJIn, dVdJOut] = ReadLevelInfo(iT, EeV);
%       for iMol = 1:NMolecules
%         EeV(:,iMol) = EeV(:,iMol) + LevelEeV(1,iMol);
%       end
%     %end
% 
%  
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% READING RATES
%     if PlotProcProbFlg == 1 || PlotKsFlg == 1 || PlotKsMoleFracsFlg == 1 || PlotMovieFlg == 1 || WriteNetworkCSVFlg == 1 || ComputeMasterEquationFlg == 1 || PlotEnergyDepletionsFlg == 1
%        [RatesMatrix, DissRates, ProcessesRates, KRec, KRecOverg] = ReadRates(iT, RatesMatrix, DissRates, ProcessesRates, StartBin, FinalBin, DissCorrectionFactor, NLevels, Levelg, LevelEeV);
%     
%        [DissRates, ProcessesRates] = CorrectRatesCG(iT, TypeCG, NLevels, Levelvqn, Levelg, LevelEeV0, DeltaEintDiss, DissRates, ProcessesRates);
%     end
%     
% %     figure(1)
% %     semilogy(LevelEeV,DissRates, 'og')
% %     hold on
%     if (DissFlg > 1)
%       
%       opts = delimitedTextImportOptions("NumVariables", 7);
%       opts.DataLines = [1, Inf];
%       opts.Delimiter = ["(", ")", ",", ":"];
%       opts.VariableNames = ["Var1", "VarName2", "Var3", "E12", "Var5", "Var6", "Var7"];
%       opts.SelectedVariableNames = ["VarName2", "E12"];
%       opts.VariableTypes = ["string", "double", "string", "double", "string", "string", "string"];
%       opts = setvaropts(opts, [1, 3, 5, 6, 7], "WhitespaceRule", "preserve");
%       opts = setvaropts(opts, [1, 3, 5, 6, 7], "EmptyFieldRule", "auto");
%       opts.ExtraColumnsRule = "ignore";
%       opts.EmptyLineRule = "read";
%       if DissFlg == 2 || DissFlg == 6
%         tbl = readtable(strcat("/home/venturi/WORKSPACE/Mars_Database/Run_0D/database/kinetics/O3_UMN/T",num2str((T0_Vec(iT))),"K/Diss_Corrected.dat"), opts);
%       elseif DissFlg == 3
%         tbl = readtable(strcat("/home/venturi/WORKSPACE/Mars_Database/Run_0D/database/kinetics/O3_UMN/T",num2str((T0_Vec(iT))),"K/Diss_VS.dat"), opts);
%       elseif DissFlg == 4
%         tbl = readtable(strcat("/home/venturi/WORKSPACE/Mars_Database/Run_0D/database/kinetics/O3_UMN/T",num2str((T0_Vec(iT))),"K/Diss_Phys_",num2str(TempBins),"Bins.dat"), opts);
%       elseif DissFlg == 5
%         tbl = readtable(strcat("/home/venturi/WORKSPACE/Mars_Database/Run_0D/database/kinetics/O3_UMN/T",num2str((T0_Vec(iT))),"K/Diss_Phys_Fitted_",num2str(TempBins),"Bins.dat"), opts);
%       end 
%       Idx      = tbl.VarName2;
%       DissTemp = tbl.E12;
%       clear opts tbl
%       ProcessesRates               = zeros(NLevels,4,1);
%       ProcessesRates(Idx(:,1),1,1) = DissTemp(:,1);
%       DissRates                    = ProcessesRates(:,1,1);
% %       semilogy(LevelEeV,DissRates, 'or')
%     end
%     
%     
%     %% READING AND COMPUTING QUANTITIES
%     [t, MolFracs, Ttrans, rho, P, nd, E] = ReadBox(iT); 
% 
% 
%     if PlotTemperaturesFlg == 1
%       [Tint] = ReadInternalT(iT);
%     end
% 
%     
%     if PlotProcProbFlg == 1 || PlotKsFlg == 1 || PlotKsMoleFracsFlg == 1 || PlotEnergiesFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || WriteNetworkCSVFlg == 1 || ComputeMasterEquationFlg == 1 || PlotEnergyDepletionsFlg == 1
%       if ReUpload == 1
%         [PopVec, PopOverQ, Pop, NSteps] = ReadPop(iT, QBins);
%         [niRatio, ProcessesRatesOverall] = ComputeProcessesRatesOverall(iT, Pop, ProcessesRates);
%         filenamePop = strcat(OutputPath,'/Pop')
%         save(filenamePop, 'PopVec', 'PopOverQ', 'Pop', 'NSteps', 'niRatio', 'ProcessesRatesOverall', '-v7.3')
%       elseif ReUpload == 2
%         clear PopVec PopOverQ Pop NSteps niRatio ProcessesRatesOverall
%         filenamePop = strcat(OutputPath,'/Pop')
%         load(filenamePop, 'PopVec', 'PopOverQ', 'Pop', 'NSteps', 'niRatio', 'ProcessesRatesOverall')
%         %FileInelDissEigs = strcat('./FileInelDissEigs')
%         %load(FileInelDissEigs, 'A', 'D', 'DInv', 'R', 'L')
%       end
%       
%       [tStart, iQSS, tEnd] =  FindQSS(iT, t, ProcessesRatesOverall);
%     end
%     
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% COMPUTING QUANTITIES
%     if PlotProcProbFlg == 1 || PlotEnergiesFlg == 1 || PlotQSSDistrFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || WriteNetworkCSVFlg == 1 || ComputeMasterEquationFlg == 1 || PlotEnergyDepletionsFlg == 1
%       [StpInstants] = ComputeStepInstants(t);        
%     end
%     
%       
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% READING / COMPUTING LEVELS-To-BINs MAPPING
%     if LevToBinFlg == 1
%       %[LevToBin, qnToBin] = ReadLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, StpInstants, GroupNb);
%       BinEnergy    = 0.d0 .* EeV;
%       BinEnergyMax = 0.d0 .* EeV;
%       opts = delimitedTextImportOptions("NumVariables", 2);
%       opts.DataLines = [2, Inf];
%       opts.Delimiter = ",";
%       opts.VariableNames = ["Var1", "Bin"];
%       opts.SelectedVariableNames = "Bin";
%       opts.VariableTypes = ["string", "double"];
%       opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
%       opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
%       opts.ExtraColumnsRule = "ignore";
%       opts.EmptyLineRule = "read";
%       tbl = readtable("/home/venturi/WORKSPACE/SpectralCluster/output/T10000K/11Bins_wExchT_InelDiss.dat", opts);
%       LevToBin = tbl.Bin;
%       clear opts tbl
%     elseif LevToBinFlg == 2
%       [LevToBin, BinEnergy, BinEnergyMax, qnToBin] = ComputeLevToBin(NLevels, Levelvqn, Leveljqn, LevelEeV, DeltaEintDiss, DeltaEintDepth);
%     else
%       LevToBin     = 0.d0 .* LevelEeV + 1;
%       BinEnergy    = 0.d0 .* EeV;
%       BinEnergyMax = 0.d0 .* EeV;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
% 
%     if PlotProcProbFlg == 1 || PlotEnergiesFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || PlotQSSDistrFlg == 1 || PlotEnergyDepletionsFlg == 1
%       if PlotMovieFlg == 1 || PlotStepsFlg == 1 || PlotEnergiesFlg == 1 || PlotEnergyDepletionsFlg == 1
%         [Steps] = ComputeMovieSteps(iT, t);  
%       end
% 
%       if PlotEnergiesFlg == 1 || PlotEnergyDepletionsFlg == 1
%         [eInt, eRot, eVib] = ComputeEnergies(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps);
%       end
%       
%       if PlotEnergyDepletionsFlg == 1
%         [CDInt, CDVib, CDRot] = ComputeEnergyDepletions(NLevels, Levelg, LevelEeV0, LevelEeV, LevelEeVVib0, LevelEeVRot, Pop, ProcessesRates, MolFracs, nd, LevelToBin, QBins, Steps);
%       end
% 
%       if PlotQSSDistrFlg == 1
%          [qssInt, qssRot, qssVib] = ComputeQssDistr(NLevels, LevelQ, Levelg, Pop, QBins, LevelEeV0, vEeVVib, vEeVVib0, LevelEeVVib0, LevelEeVRot, Levelvqn, LevelToBin, Steps); 
%       end
%     else
%       vEeVVib = [];
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     
%     ReUpload = 0
%   end
%    
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %% PLOTTING
% 
%   % Plotting MOLE FRACTIONS
%   if PlotMolFracsFlg == 1
%     [iFigure] = PlotMolFracs(iT, iFigure, t, MolFracs);     
%   end
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
% 
%   % Plotting TEMPERATURES
%   if PlotTemperaturesFlg == 1
%     [iFigure] = PlotTemperatures(iT, iFigure, t, Tint, Ttrans);
%   end 
% 
% 
%   % Plotting QSS RATES
%   if PlotKsFlg == 1 
%     [iFigure] = PlotKs(iT, iFigure, t, ProcessesRatesOverall);  
%   end
% 
% 
%   % QSS RATES PLUS MOLE FRACTION
%   if PlotKsMoleFracsFlg == 1 
%     [iFigure] = PlotKsMoleFracs(iT, iFigure, t, ProcessesRatesOverall, MolFracs, tStart, iQSS, tEnd);                    
%   end
% 
%   
%   % Plotting MOVIES
%   if PlotMovieFlg == 1
%     iFigure = PlotMovie(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, DeltaEintDiss, LevelQ, Levelg, Levelvqn, Leveljqn, LevelToBin, Pop, QBins, Steps, LevToBin, ProcessesRates);  
%   end
%   
% 
%   % Plotting STEPS
%   if PlotStepsFlg == 1
%     [iFigure] = PlotSteps(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, LevelEeV0, LevelQ, Levelg, LevelToBin, Pop, QBins, StpInstants, Levelvqn, Leveljqn, LevToBin, DeltaEintDiss, ProcessesRates);
%     %[iFigure] = PlotPartFunc(iT, iFigure, t, MolFracs, ProcessesRatesOverall, NLevels, LevelEeV, LevelEeV0, LevelQ, Levelg, LevelToBin, Pop, QBins, StpInstants, Levelvqn, Leveljqn, LevToBin, DeltaEintDiss, ProcessesRates)    
%   end
%   
%   
%   % Plotting Energies
%   if PlotEnergiesFlg == 1
%     [iFigure] = PlotEnergies(iT, iFigure, t, eInt, eRot, eVib, MolFracs, P, Steps);
%   end
%   
%   
%   % Plotting Energy Deplations
%   if PlotEnergyDepletionsFlg == 1
%     [iFigure] = PlotEnergyDepletions(iFigure, t, CDInt, CDVib, CDRot, tStart, iQSS, tEnd);
%   end
%   
%   
%   % Plotting Processes Probabilities
%   if PlotProcProbFlg == 1
%     [iFigure] = PlotProcProb(iT, iFigure, t, ProcessesRates, ProcessesRatesOverall, niRatio, LevelEeV, DeltaEintDiss, vEeVVib, Steps, StpInstants, Leveljqn, Levelvqn);
%   end
%   
%   
%   % Reconstructing Q.B. Levels from updated Saha Eq.
%   if PlotNewSahaFlg == 1
%     [iFigure] = PlotNewSaha(iT, iFigure, iFigureStart, StpInstants, NLevels, Levelg, LevelEeV, DeltaEintDiss, LevToBin, BinEnergy, BinEnergyMax, t, MolFracs, nd, VMax, ProcessesRates, RatesMatrix, Tint)
%   end 
%   
%   
%   if PlotTausFlg == 1
%     [iFigure] = PlotTaus(iFigure);
%   end
%   
%   
%   if WriteNetworkCSVFlg == 1
%     WriteNetworkCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, Pop, StpInstants, LevelEeVVib0, LevelEeVRot)
%   end
%   
%   
%   if WriteRatesToParaviewCSVFlg == 1
%     WriteRatesToParaviewCSV(iT, NLevels, LevelEeV, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, LevelEeVVib0, LevelEeVRot)
%   end
%   
%   if ComputeMasterEquationFlg == 1
%    ComputeMasterEquation(iT, NLevels, LevelEeV, LevelEeV0, Levelvqn, Leveljqn, rIn, rOut, Levelg, RatesMatrix, ProcessesRates, t, MolFracs, nd, LevelEeVVib0, LevelEeVRot)
%   end
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   
% end