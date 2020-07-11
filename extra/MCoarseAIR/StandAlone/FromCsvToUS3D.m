close all
clear all
clc


NBins1    = 6115;
NBins2    = 6115;
TVec      = [1500.0, 2500.0, 5000.0, 6000.0, 8000.0, 10000.0, 12000.0, 12000.0, 14000.0, 15000.0, 20000.0];
NT        = length(TVec);
% RatesFldr = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/kinetics/O3_UMN_DP10/';
% DestFldr  = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/O3_UMN/';
RatesFldr = '/home/venturi/Desktop/O3_UMN/O3_UMN/';
DestFldr  = '/home/venturi/Desktop/O3_UMN/';
HomoVec   = [1];
HeteroVec = [];

% NBins1    = 10;
% NBins2    = 10;
% TVec      = [2500.0, 3500.0, 7500.0, 10000.0, 12500.0, 15000.0, 20000.0, 25000.0, 30000.0, 40000.0, 50000.0];
% NT        = length(TVec);
% RatesFldr = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/kinetics/N3_NASA_DP10/';
% DestFldr  = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/N3_NASA/';
% HomoVec   = [];
% HeteroVec = [];

% NBins1    = 10;
% NBins2    = 10;
% TVec      = [2500.0, 5000.0, 7500.0, 10000.0, 15000.0, 20000.0];
% NT        = length(TVec);
% RatesFldr = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/kinetics/N2O_UMN_DP10/';
% DestFldr  = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/N2O_UMN/';
% HomoVec   = [];
% HeteroVec = [1];

% NBins1    = 10;
% NBins2    = 10;
% TVec      = [2500.0, 5000.0, 7500.0, 10000.0, 15000.0, 20000.0];
% NT        = length(TVec);
% RatesFldr = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/kinetics/NON_UMN_DP10/';
% DestFldr  = '/home/venturi/WORKSPACE/Air_Database/US3D_Database/NON_UMN/';
% HomoVec   = [2];
% HeteroVec = [1];


FormatT = '%e';
TVecChar = num2str(TVec(1));
for iT=2:NT
   FormatT  = strcat(FormatT, ' %e');
   TVecChar = strcat(TVecChar, ',', num2str(TVec(iT)));
end

iT = 1;
KDMat     = zeros(NBins1, NT);
KExcitMat = zeros(NBins1, NBins1, NT);
KExchMat  = zeros(NBins1, NBins2, NT);
for TT = TVec
    RatesFldrT = strcat(RatesFldr, '/T', num2str(floor(TT)) , 'K/');
    
    
    opts = delimitedTextImportOptions("NumVariables", 2);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["P", "KD"];
    opts.VariableTypes = ["double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    tbl = readtable(strcat(RatesFldrT,'/Diss_Corrected.csv'), opts);
    %tbl = readtable(strcat(RatesFldrT,'/Diss.csv'), opts);
    P  = tbl.P;
    KD = tbl.KD;
    clear opts tbl
    KDMat(P(:),iT) = KD(:);
    clear P KD
 
    
%     opts = delimitedTextImportOptions("NumVariables", 3);
%     opts.DataLines = [2, Inf];
%     opts.Delimiter = ",";
%     opts.VariableNames = ["P1", "Q", "k_P"];
%     opts.VariableTypes = ["double", "double", "double"];
%     opts.ExtraColumnsRule = "ignore";
%     opts.EmptyLineRule = "read";
%     tbl = readtable(strcat(RatesFldrT,'/Inel.csv'), opts);
%     P     = tbl.P1;
%     Q     = tbl.Q;
%     KInel = tbl.k_P;
%     clear opts tbl
%     for ii=1:length(P)
%         KExcitMat(P(ii),Q(ii),iT) = KInel(ii);
%     end
%     clear P Q kInel
    
%     
%     for iHomo = HomoVec
% 
%         opts = delimitedTextImportOptions("NumVariables", 3);
%         opts.DataLines = [2, Inf];
%         opts.Delimiter = ",";
%         opts.VariableNames = ["P1", "Q", "k_P"];
%         opts.VariableTypes = ["double", "double", "double"];
%         opts.ExtraColumnsRule = "ignore";
%         opts.EmptyLineRule = "read";
%         tbl = readtable(strcat(RatesFldrT,'/Exch_Type',num2str(iHomo),'.csv'), opts);
%         P     = tbl.P1;
%         Q     = tbl.Q;
%         KHomo = tbl.k_P;
%         clear opts tbl
%         for ii=1:length(P)
%             KExcitMat(P(ii),Q(ii),iT) = KExcitMat(P(ii),Q(ii),iT) + KHomo(ii);
%         end
%         clear P Q KHomo
% 
%     end
%     
%     
%     for iHetero = HeteroVec
% 
%         opts = delimitedTextImportOptions("NumVariables", 3);
%         opts.DataLines = [2, Inf];
%         opts.Delimiter = ",";
%         opts.VariableNames = ["P1", "Q", "k_P"];
%         opts.VariableTypes = ["double", "double", "double"];
%         opts.ExtraColumnsRule = "ignore";
%         opts.EmptyLineRule = "read";
%         tbl = readtable(strcat(RatesFldrT,'/Exch_Type',num2str(iHetero),'.csv'), opts);
%         P       = tbl.P1;
%         Q       = tbl.Q;
%         KHetero = tbl.k_P;
%         clear opts tbl
%         for ii=1:length(P)
%             KExchMat(P(ii),Q(ii),iT) = KExchMat(P(ii),Q(ii),iT) + KHetero(ii);
%         end
%         clear P Q KHetero
% 
%     end
    
    
    iT = iT + 1;
end



FileName1 = strcat(DestFldr, '/Dissociation.csv');
fileID1   = fopen(FileName1,'w');
fprintf(fileID1,strcat('#T=',TVecChar,'\n'));
for i=1:NBins1
   if (sum(KDMat(i,:)) > 0.0)
        fprintf(fileID1,['%i ',FormatT,'\n'], i, KDMat(i,:));
   end
end
fclose(fileID1);


% FileName1 = strcat(DestFldr, '/Excitation.csv');
% fileID1   = fopen(FileName1,'w');
% fprintf(fileID1,strcat('#T=',TVecChar,'\n'));
% for i=1:NBins1
%    for j=1:NBins1
%        if (sum(KExcitMat(i,j,:)) > 0.0)
%             fprintf(fileID1,['%i %i ',FormatT,'\n'], i, j, KExcitMat(i,j,:));
%        end
%    end
% end
% fclose(fileID1);
% 
% 
% FileName1 = strcat(DestFldr, '/Exchange.csv');
% fileID1   = fopen(FileName1,'w');
% fprintf(fileID1,strcat('#T=',TVecChar,'\n'));
% for i=1:NBins1
%    for j=1:NBins2
%        if (sum(KExchMat(i,j,:)) > 0.0)
%             fprintf(fileID1,['%i %i ',FormatT,'\n'] , i, j, KExchMat(i,j,:));
%        end
%    end
% end
% fclose(fileID1);
