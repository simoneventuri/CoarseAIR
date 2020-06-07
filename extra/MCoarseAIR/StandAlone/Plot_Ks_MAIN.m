close all
clear all
clc


global Param Input Syst
iFig = 1;

Input.SystNameLong       = 'CO2_NASA'
Input.Paths.MainFldr     = '/home/venturi/WORKSPACE/Mars_Paper/'
Input.Paths.KGlobal      = strcat(Input.Paths.MainFldr, '/Data/', Input.SystNameLong, '/KGlobal_9_1_1_1_O2.csv')
NPESs                    = 0
ExchToMol                = [2];

Input.FigureFormat       = 'PrePrint';
Input.iFig               = 101;
Input.SaveFigsFlgInt     = 0;
Input.Paths.SaveFigsFldr = strcat(Input.Paths.MainFldr, '/Figures/', Input.SystNameLong);


Syst.NameLong = Input.SystNameLong;
Syst          = Initialize_ChemicalSyst(Syst)
Initialize_Parameters()


Plot_Ks(ExchToMol)

% figure(4321)
% for iPES = 1:NPESs
%     
%     opts.DataLines = [2, Inf];
%     opts.Delimiter = ",";
%     if (Syst.NProc == 3)
%         opts = delimitedTextImportOptions("NumVariables", 3);
%         opts.VariableNames = ["TVec", "KDEq", "KExEq"];
%         opts.VariableTypes = ["double", "double", "double"];
%     else
%         opts = delimitedTextImportOptions("NumVariables", 4);
%         opts.VariableNames = ["TVec", "KDEq", "KExEq1", "KExEq"];
%         opts.VariableTypes = ["double", "double", "double", "double"];
%     end
%     opts.ExtraColumnsRule = "ignore";
%     opts.EmptyLineRule = "read";
%     FileName = strcat(Input.Paths.MainFldr, "/Data/", Input.SystNameLong, "_PES", num2str(iPES), "/KEq.csv");
%     tbl = readtable(FileName, opts);
%     TVec  = tbl.TVec;
%     KDEq  = tbl.KDEq;
%     KExEq = tbl.KExEq;
%     clear opts tbl
%     KTotEq = KExEq + KDEq; 
%     
%     semilogy(10000./TVec, KDEq,   '.')
%     hold on
%     semilogy(10000./TVec, KExEq,  's')
% %     semilogy(10000./TVec, KTotEq, 'x')
% 
% end