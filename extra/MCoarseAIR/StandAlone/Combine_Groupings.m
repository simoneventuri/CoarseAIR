close all
clear all
clc


% NLevels  = 6115;
% InelFile = "/home/venturi/WORKSPACE/SpectralCluster/data/O3_UMN/T10000K/GroupingInfo_6Bins_WExch.csv";
% CBFile   = "/home/venturi/WORKSPACE/SpectralCluster/data/O3_UMN/T10000K/LevelsInfo_CB.csv";
% DestFldr = "/home/venturi/WORKSPACE/Air_Database/Run_0D/database/grouping/";

NLevels  = 9390;
InelFile = "/home/venturi/WORKSPACE/SpectralCluster/data/N3_NASA/T10000K/GroupingInfo_30Bins_WExch.csv";
CBFile   = "/home/venturi/WORKSPACE/SpectralCluster/data/N3_NASA/T10000K/LevelsInfo_CB.csv";
DestFldr = "/home/venturi/WORKSPACE/Air_Database/Run_0D/database/grouping/";





Group    = zeros(NLevels,1);


opts = delimitedTextImportOptions("NumVariables", 7);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "iLevel", "Group"];
opts.SelectedVariableNames = ["iLevel", "Group"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5"], "EmptyFieldRule", "auto");
tbl = readtable(InelFile, opts);
iLevel_Inel = tbl.iLevel;
Group_Inel  = tbl.Group;
clear opts tbl
NbGroups_Inel      = max(Group_Inel);
Group(iLevel_Inel) = Group_Inel;


opts = delimitedTextImportOptions("NumVariables", 9);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Group_CB", "iLevel"];
opts.SelectedVariableNames = ["Group_CB", "iLevel"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7"], "EmptyFieldRule", "auto");
tbl = readtable(CBFile, opts);
Group_CB  = tbl.Group_CB;
iLevel_CB = tbl.iLevel;
clear opts tbl
NbGroups_CB      = max(Group_CB);
Group(iLevel_CB) = Group_CB + NbGroups_Inel;


FileName1 = strcat(DestFldr,"/", num2str(NbGroups_Inel), "Inel_", num2str(NbGroups_CB), "CB.csv");
fileID1   = fopen(FileName1,'w');
fprintf(fileID1,'#Idx,Group\n');
for iLevel=1:NLevels
    fprintf(fileID1,'%i,%i\n', iLevel, Group(iLevel));
end
fclose(fileID1);