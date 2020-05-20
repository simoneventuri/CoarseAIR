close all
clear all
clc

CoarseAIR_RatesFldr = '/home/venturi/WORKSPACE/CoarseAIR/O2C_ALL_OLD/Test/CO2/Rates_TEMP/';
CGQCT_RatesFldr     = '/home/venturi/WORKSPACE/CG-QCT/run_O2C_ALL/Test/CO2/O2/Rates/';
TTran               = 12500;

NLevels1   = 6078;
NLevels2   = 13521;
NLevels3   = 13521;


FileWithProcesses = "/home/venturi/WORKSPACE/CoarseAIR/TTTTEMP/O2C_ALL_TEMP/input/ProcessesToRunList.inp";
opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = "VarName1";
opts.VariableTypes = "double";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
ProcessesToRunList = readtable(FileWithProcesses, opts);
clear opts
iLevelVec = ProcessesToRunList.VarName1;


for iLevel = iLevelVec'
    iLevel
    
    RatesFile = strcat(CoarseAIR_RatesFldr, "/T_", num2str(TTran), "_", num2str(TTran), "/Proc", num2str(iLevel), ".csv");
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [6, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Process", "Rate", "RateSD"];
    opts.VariableTypes = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(RatesFile, opts);
    Process = tbl.Process;
    Rate    = tbl.Rate;
    RateSD  = tbl.RateSD;
    clear opts tbl
    CoarseAIR_Rates   = zeros(NLevels1+NLevels2+NLevels3+4);
    CoarseAIR_RatesSD = zeros(NLevels1+NLevels2+NLevels3+4);
    
    
    CoarseAIR_Rates(Process(:))   = Rate(:);
    CoarseAIR_RatesSD(Process(:)) = RateSD(:);
    
    
    CGQCT_Rates   = zeros(NLevels1+NLevels2+NLevels3+1);
    CGQCT_RatesSD = zeros(NLevels1+NLevels2+NLevels3+1);
    
    CGQCT_Rates(1)                                                  = CoarseAIR_Rates(2) + CoarseAIR_Rates(NLevels1+3) + CoarseAIR_Rates(NLevels1+NLevels2+4);
    CGQCT_Rates(2:NLevels1+1)                                       = CoarseAIR_Rates(3:NLevels1+2);
    CGQCT_Rates(NLevels1+2:NLevels1+NLevels2+1)                     = CoarseAIR_Rates(NLevels1+4:NLevels1+NLevels2+3);
    CGQCT_Rates(NLevels1+NLevels2+2:NLevels1+NLevels2+NLevels3+1)   = CoarseAIR_Rates(NLevels1+NLevels2+5:NLevels1+NLevels2+NLevels3+4);
    
    CGQCT_RatesSD(1)                                                = sqrt( CoarseAIR_RatesSD(2)^2 + CoarseAIR_RatesSD(NLevels1+3)^2 + CoarseAIR_RatesSD(NLevels1+NLevels2+4)^2 );
    CGQCT_RatesSD(2:NLevels1+1)                                     = CoarseAIR_RatesSD(3:NLevels1+2);
    CGQCT_RatesSD(NLevels1+2:NLevels1+NLevels2+1)                   = CoarseAIR_RatesSD(NLevels1+4:NLevels1+NLevels2+3);
    CGQCT_RatesSD(NLevels1+NLevels2+2:NLevels1+NLevels2+NLevels3+1) = CoarseAIR_RatesSD(NLevels1+NLevels2+5:NLevels1+NLevels2+NLevels3+4);
    
    
    RatesFile = strcat(CGQCT_RatesFldr, "/T_", num2str(TTran), "_", num2str(TTran), "/Bin", num2str(iLevel), ".dat");
    fileID    = fopen(RatesFile,'w');
    fprintf(fileID,'#\n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'#\n');
    fprintf(fileID,'#\n');

    for iProc = 1:NLevels1+NLevels2+NLevels3+1
        if (CGQCT_Rates(iProc) > 0.0)
            fprintf(fileID,'--%38i%20.10e%20.10e\n', iProc, CGQCT_Rates(iProc), CGQCT_RatesSD(iProc));
        end
    end
    fclose(fileID);
    
    
    clear Process Rate RateSD CoarseAIR_Rates CoarseAIR_RatesSD CGQCT_Rates CGQCT_RatesSD
end 