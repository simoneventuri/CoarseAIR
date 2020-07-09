close all
clear all
clc

TestFldr = "/home/venturi/WORKSPACE/CoarseAIR/O4_UMN_CGM_10Groups/Test/";
%TestFldr = "/home/venturi/WORKSPACE/CoarseAIR/O4_TEST/Test/";
Vel      = 0.1018961233D-10;



FileName = strcat(TestFldr, "/OaObOcOd/OaOb/Bins_10/QNsEnBin.csv");
opts = delimitedTextImportOptions("NumVariables", 6);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "Levelvqn", "Leveljqn", "LevelEEh", "Levelg", "LevelToBin"];
opts.SelectedVariableNames = ["Levelvqn", "Leveljqn", "LevelEEh", "Levelg", "LevelToBin"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
tbl = readtable(FileName, opts);
Levelvqn   = tbl.Levelvqn;
Leveljqn   = tbl.Leveljqn;
LevelEEh   = tbl.LevelEEh;
Levelg     = tbl.Levelg;
LevelToBin = tbl.LevelToBin;
clear opts tbl
NLevels = length(LevelToBin);
Nvqn    = max(Levelvqn)+1; 
Njqn    = max(Leveljqn)+1;

%LevelToBin = [1:1:NLevels];
NBins   = max(LevelToBin);

QNToBin = ones(Njqn+20,Nvqn+20).*-1.0;
for iLevel=1:NLevels
    QNToBin(Leveljqn(iLevel)+1,Levelvqn(iLevel)+1) = LevelToBin(iLevel);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking Statistics
%%% 
% NArr = 64;
% NTot = Nvqn^2*Njqn^2*NArr;
% ii=1;
% for kBin=10:10%NBins
%     for lBin=kBin:10%NBins
%         fprintf('iBin = %i, jBin = %i \n', kBin, lBin);
%         
% 
%         FileName = strcat(TestFldr, "/T_10000_10000/Bins_", num2str(kBin), "_", num2str(lBin), "/trajectories.csv");
%         opts = delimitedTextImportOptions("NumVariables", 14);
%         opts.DataLines = [2, Inf];
%         opts.Delimiter = ",";
%         opts.VariableNames = ["Var1", "Var2", "bmax", "b_i", "Var5", "Var6", "Var7", "Var8", "Var9", "j1_f", "v1_f", "j2_f", "v2_f", "arr_f"];
%         opts.SelectedVariableNames = ["bmax", "b_i", "j1_f", "v1_f", "j2_f", "v2_f", "arr_f"];
%         opts.VariableTypes = ["string", "string", "double", "double", "string", "string", "string", "string", "string", "double", "double", "double", "double", "double"];
%         opts.ExtraColumnsRule = "ignore";
%         opts.EmptyLineRule = "read";
%         opts = setvaropts(opts, ["Var1", "Var2", "Var5", "Var6", "Var7", "Var8", "Var9"], "WhitespaceRule", "preserve");
%         opts = setvaropts(opts, ["Var1", "Var2", "Var5", "Var6", "Var7", "Var8", "Var9"], "EmptyFieldRule", "auto");
%         tbl = readtable(FileName, opts);
%         bMaxTemp = tbl.bmax;
%         b_i      = tbl.b_i;
%         j1_f  = round(tbl.j1_f  - 0.5);
%         v1_f  = round(tbl.v1_f  - 0.5);
%         j2_f  = round(tbl.j2_f  - 0.5);
%         v2_f  = round(tbl.v2_f  - 0.5);
%         arr_f = round(tbl.arr_f - 0.5);
%         clear opts tbl
%         bMax  = unique(bMaxTemp);
%         Nbb   = length(bMax);
%         
%         IdVec   = [];
%         NRepeat = [];
%         for ib=1:Nbb
%             b(ib).IdVec    = [];
%         end
%         ARing(1) = bMax(1)^2 * pi;
%         for ib=2:length(bMax)
%             ARing(ib) = (bMax(ib)^2 - bMax(ib-1)^2) * pi;
%         end
% 
%         bVec     = zeros(length(bMax),1);
%         for iTraj = 1:length(arr_f)
%             if (arr_f(iTraj) > 0)
% 
%                 ib       = sum(b_i(iTraj) > bMax) + 1;
%                 bVec(ib) = bVec(ib) + 1;
% 
%                 iP     = floor(arr_f(iTraj)/16);
%                 iType  =   mod(arr_f(iTraj),16);
%                 iTypeA = floor(   mod(iType,4) / 2.0);
%                 iTypeB = floor( floor(iType/4) / 2.0);
% 
%                 iBin = 0;
%                 if (iTypeA == 0)
%                     iBin = QNToBin(j1_f(iTraj)+1,v1_f(iTraj)+1); 
%                     iBin = iBin + 1;
%                 end
%                 jBin = 0;
%                 if (iTypeB == 0)
%                     jBin = QNToBin(j2_f(iTraj)+1,v2_f(iTraj)+1); 
%                     jBin = jBin + 1;
%                 end
% 
%                 Id   =      v1_f(iTraj)*Njqn*Nvqn*Njqn*NArr;
%                 Id   = Id + j1_f(iTraj)*Nvqn*Njqn*NArr;
%                 Id   = Id + v2_f(iTraj)*Njqn*NArr;
%                 Id   = Id + j2_f(iTraj)*NArr;
%                 Id   = Id + arr_f(iTraj);
%                     
%                 jProc = 0;
%                 iProc = 1;
%                 while (jProc==0) && (iProc <= length(IdVec)) 
%                     if (Id == IdVec(iProc))
%                         jProc = iProc;
%                     end
%                     iProc=iProc+1;
%                 end
%                 if jProc > 0
%                     NRepeat(jProc,ib) = NRepeat(jProc,ib) + 1;
%                 else
%                     IdVec           = [IdVec;  Id];
%                     NRepeat         = [NRepeat; zeros(1,Nbb)];
%                     NRepeat(end,ib) = 1;
%                 end
%                 
%             end
%         end
%         
%         CrossVec = zeros(length(IdVec),1);
%         for ib=1:length(bMax)
%             b(ib).Cross(:) = NRepeat(:,ib) ./ bVec(ib) .* ARing(ib);
%             CrossVec       = CrossVec  + b(ib).Cross(:); 
%         end
%         
%         [IdSort, Order] = sort(IdVec);
%         CrossVecSort    = CrossVec(Order);
%         NewCross        = [IdSort, CrossVecSort];
%         
%         
%     
%         FileName = strcat(TestFldr, "/T_10000_10000/Bins_", num2str(kBin), "_", num2str(lBin), "/statistics.csv");
%         startRow = 2;
%         formatSpec = '%6f%8f%8f%8f%8f%*8*s%*8*s%*8*s%*8*s%*8f%18f%[^\n\r]';
%         fileID = fopen(FileName,'r');
%         dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%         fclose(fileID);
%         j1_ff    = dataArray{:, 1};
%         v1_ff    = dataArray{:, 2};
%         j2_ff    = dataArray{:, 3};
%         v2_ff    = dataArray{:, 4};
%         arr_ff   = dataArray{:, 5};
%         Cross_ff = dataArray{:, 6};
%         clearvars filename startRow formatSpec fileID dataArray ans;
%         IdVecOld    = [];
%         CrossVecOld = [];
%         AA1=Njqn*Nvqn*Njqn*NArr;
%         AA2=Nvqn*Njqn*NArr;
%         AA3=Njqn*NArr;
%         AA4=NArr;
%         for i=1:length(j1_ff)
%             Id          = v1_ff(i)*AA1 + j1_ff(i)*AA2 + v2_ff(i)*AA3 + j2_ff(i)*AA4 + arr_ff(i);            
%             IdVecOld    = [IdVecOld; Id];
%             CrossVecOld = [CrossVecOld; Cross_ff(i)];
%         end
%         [IdSortOld, Order] = sort(IdVecOld);
%         CrossVecSortOld    = CrossVecOld(Order);
%         OldCross           = [IdSortOld, CrossVecSortOld];
%         
%         
%         for i=1:size(OldCross,1)
%             Error1 = abs(OldCross(i,1)-NewCross(i,1));
%             Error2 = abs(OldCross(i,2)-NewCross(i,2));
%             if Error1 > 0 || Error2 > 1e-5
%                 clear v1_f j1_f v2_f j2_f arr_f v1_ff j1_ff v2_ff j2_ff arr_ff
%                 
%                 v1_f   = floor(NewCross(i,1)/ (AA1));
%                 Tempp  =   mod(NewCross(i,1), AA1);
%                 j1_f   = floor(Tempp/ (AA2));
%                 Tempp  =   mod(Tempp, AA2);
%                 v2_f   = floor(Tempp/ (AA3));
%                 Tempp  =   mod(Tempp, AA3);
%                 j2_f   = floor(Tempp/ AA4);
%                 Tempp  =   mod(Tempp, AA4);
%                 arr_f  = Tempp;
%                 fprintf('vqn1 = %i; jqn1 = %i; vqn2 = %i; jqn2 = %i; arr_f = %i; Cross = %e\n', v1_f, j1_f, v2_f, j2_f, arr_f, NewCross(i,2))
% 
%                 v1_ff  = floor(OldCross(i,1)/ (Njqn*Nvqn*Njqn*NArr));
%                 Tempp  =   mod(OldCross(i,1), Njqn*Nvqn*Njqn*NArr);
%                 j1_ff  = floor(Tempp/ (Nvqn*Njqn*NArr));
%                 Tempp  =   mod(Tempp, Nvqn*Njqn*NArr);
%                 v2_ff  = floor(Tempp/ (Njqn*NArr));
%                 Tempp  =   mod(Tempp, Njqn*NArr);
%                 j2_ff  = floor(Tempp/ NArr);
%                 Tempp  =   mod(Tempp, NArr);
%                 arr_ff = Tempp;
%                 fprintf('vqn1 = %i; jqn1 = %i; vqn2 = %i; jqn2 = %i; arr_f = %i; Cross = %e\n\n', v1_ff, j1_ff, v2_ff, j2_ff, arr_ff, OldCross(i,2))
%                 
%                 pause
%             end
%         end
%         
%         
%         ii=ii+1;
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking Postprocessing
%%% 
ErrorVec = zeros(NBins*(NBins+1)/2,1);
ii=1;
for kBin=1:NBins
    for lBin=kBin:NBins
        fprintf('iBin = %i, jBin = %i \n', kBin, lBin);
        

        FileName = strcat(TestFldr, "/OaObOcOd/Rates/T_10000_10000/i", num2str(kBin), "_j", num2str(lBin), ".csv");
        opts = delimitedTextImportOptions("NumVariables", 3);
        opts.DataLines = [6, Inf];
        opts.Delimiter = ",";
        opts.VariableNames = ["Idx", "VarName2", "Var3"];
        opts.SelectedVariableNames = ["Idx", "VarName2"];
        opts.VariableTypes = ["double", "double", "string"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        opts = setvaropts(opts, "Var3", "WhitespaceRule", "preserve");
        opts = setvaropts(opts, "Var3", "EmptyFieldRule", "auto");
        tbl = readtable(FileName, opts);
        Idx       = tbl.Idx;
        RatesTemp = tbl.VarName2;
        clear opts tbl
        RatesOrig         = zeros(1 + 3*(NBins+1)*(NBins+1),1);
        RatesOrig(Idx(:)) = RatesTemp(:); 



        FileName = strcat(TestFldr, "/T_10000_10000/Bins_", num2str(kBin), "_", num2str(lBin), "/statistics.csv");
        startRow = 2;
        formatSpec = '%6f%8f%8f%8f%8f%*8*s%*8*s%*8*s%*8*s%*8f%18f%[^\n\r]';
        fileID = fopen(FileName,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        j1_ff    = dataArray{:, 1};
        v1_ff    = dataArray{:, 2};
        j2_ff    = dataArray{:, 3};
        v2_ff    = dataArray{:, 4};
        arr_ff   = dataArray{:, 5};
        Cross_ff = dataArray{:, 6};
        clearvars filename startRow formatSpec fileID dataArray ans;        
        
        InelMat  = zeros(NBins+1,NBins+1);
        ExchMat1 = zeros(NBins+1,NBins+1);
        ExchMat2 = zeros(NBins+1,NBins+1);
        RatesNew = zeros(1+(NBins+1)*(NBins+1)*3,1);
        
        for iTraj = 1:length(Cross_ff)


            iP     = floor(arr_ff(iTraj)/16);
            iType  =   mod(arr_ff(iTraj),16);
            iTypeA = floor(   mod(iType,4) / 2.0);
            iTypeB = floor( floor(iType/4) / 2.0);

            
            if (iTypeA == 0)
                iBin = QNToBin(j1_ff(iTraj)+1,v1_ff(iTraj)+1); 
                if iBin < 1
                    fprintf('vqn1 = %i; jqn1 = %i; vqn2 = %i; jqn2 = %i; arr_f = %i\n', v1_ff(iTraj), j1_ff(iTraj), v2_ff(iTraj), j2_ff(iTraj), arr_ff(iTraj))
                end
            else
                iBin = 0;
            end
            jBin = 0;
            if (iTypeB == 0)
                jBin = QNToBin(j2_ff(iTraj)+1,v2_ff(iTraj)+1); 
                if jBin < 1
                    fprintf('vqn1 = %i; jqn1 = %i; vqn2 = %i; jqn2 = %i; arr_f = %i\n', v1_ff(iTraj), j1_ff(iTraj), v2_ff(iTraj), j2_ff(iTraj), arr_ff(iTraj))
                end
            else
                jBin = 0;
            end
%             Curr = 1 + (iP-1)*(NBins+1)*(NBins+1) + iBin*(NBins+1) + jBin+1;
%             if iBin > -1 && jBin > -1
%                 RatesNew(Curr) = RatesNew(Curr) + Cross_ff(iTraj)*Vel;
%             end
            %fprintf('vqn1 = %i; jqn1 = %i; iTypeA = %i; iBin = %i\n', v1_ff(iTraj), j1_ff(iTraj), iTypeA, iBin)
            %fprintf('vqn2 = %i; jqn2 = %i; iTypeB = %i; jBin = %i\n', v2_ff(iTraj), j2_ff(iTraj), iTypeB, jBin)
            %fprintf('Curr = %i\n\n', Curr)
            %pause

%             if iBin > 0 && jBin > 0
% 
%             elseif iBin > 0
%                 jBin   = 0;
%                 iTypeB = 1;
%             elseif jBin > 0
%                 iBin   = 0;
%                 iTypeA = 1;
%             else     
%                 iBin   = 0;
%                 iTypeA = 1;
%                 jBin   = 0;
%                 iTypeB = 1;
%             end
%             if iP==1
%                 InelMat(iBin+iTypeA,jBin+iTypeB)  = InelMat(iBin+iTypeA,jBin+iTypeB)  + Cross_ff(iTraj);
%             elseif iP==2
%                 ExchMat1(iBin+iTypeA,jBin+iTypeB) = ExchMat1(iBin+iTypeA,jBin+iTypeB) + Cross_ff(iTraj);
%             elseif iP==3
%                 ExchMat2(iBin+iTypeA,jBin+iTypeB) = ExchMat2(iBin+iTypeA,jBin+iTypeB) + Cross_ff(iTraj);
%             end


            if iBin > -1 && jBin > -1
                if iP==1
                    InelMat(jBin+1,iBin+1)  = InelMat(jBin+1,iBin+1)  + Cross_ff(iTraj).*Vel;
                elseif iP==2
                    ExchMat1(jBin+1,iBin+1) = ExchMat1(jBin+1,iBin+1) + Cross_ff(iTraj).*Vel;
                elseif iP==3
                    ExchMat2(jBin+1,iBin+1) = ExchMat2(jBin+1,iBin+1) + Cross_ff(iTraj).*Vel;
                end
            end
               

        end
        RatesNew     = [0.0; InelMat(:); ExchMat1(:); ExchMat2(:)];
        ErrorVec(ii) = max(max(abs(RatesOrig - RatesNew)));
        
        figure(1)
        semilogy(RatesOrig)
        hold on
        semilogy(RatesNew)
        hold off
        pause
        
        ii=ii+1;
    end
end



for iBin=2:NBins+1
    jBin=iBin;
    KK = InelMat(iBin,jBin)+InelMat(jBin,iBin);
    if (KK>0)
        fprintf('i=%i,j=%i; K=%e\n', iBin-1, jBin-1, KK)
    end
    for jBin=iBin+1:NBins+1
        KK = InelMat(iBin,jBin)+InelMat(jBin,iBin);
        if (KK>0)
            fprintf('i=%i,j=%i; K=%e\n', iBin-1, jBin-1, KK)
        end
    end
end