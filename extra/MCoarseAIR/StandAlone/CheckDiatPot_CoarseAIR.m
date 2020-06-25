close all
clear all
clc

RunFldr      = '/home/venturi/WORKSPACE/CoarseAIR/O4_TEST/'
MoleculesVec = {'OaOb'}%[{'O2'}, {'NO'}]

iMol = 1;
for Molecule = MoleculesVec

    FileName = strcat(RunFldr, "/Test/PlotPES/VDiat_From_VDiat_", Molecule, ".csv")
    opts = delimitedTextImportOptions("NumVariables", 2);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["r1", "V"];
    opts.VariableTypes = ["double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(FileName, opts);
    r1  = tbl.r1;
    V   = tbl.V;
    clear opts tbl
    
    FileName = strcat(RunFldr, "/Test/PlotPES/dVDiat_From_VDiat_", Molecule, ".csv")
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["r1", "V", "dV"];
    opts.VariableTypes = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(FileName, opts);
    r1_  = tbl.r1;
    V_   = tbl.V;
    dV_  = tbl.dV;
    clear opts tbl

    figure(iMol)
    plot(r1,V,   '-k')
    hold on
    plot(r1_,V_,  ':r')
    plot(r1_,dV_, '-r')
    
    iMol = iMol+1;
end


for iP = 1:6
   
    FileName = strcat(RunFldr, "/Test/PlotPES/dVDiat_From_PES1.csv.", num2str(iP))
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["r1", "V", "dV"];
    opts.VariableTypes = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(FileName, opts);
    r1__  = tbl.r1;
    V__   = tbl.V;
    dV__  = tbl.dV;
    clear opts tbl
    
    
    figure(10+iP)
    plot(r1__,V__,  'g')
    hold on
    plot(r1__,dV__, 'b')
    
end




FileName = strcat(RunFldr, "/Test/PlotPES/VDiat_From_VDiat_", Molecule, ".csv")
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["r1", "V"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
tbl = readtable(FileName, opts);
r1  = tbl.r1;
V   = tbl.V;
clear opts tbl
    
figure(20+1)
plot(r1,V,   '-k')
hold on
plot(r1_,V_,  ':r')
    


