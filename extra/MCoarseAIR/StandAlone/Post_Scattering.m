close all
clear all
clc


NTrajs =100000;
Dist   = 33.7e-2 * 1.889725989e+10;

%%% Physics Parameters 
AvN           = 6.0221409e+23;
AMUToKg       = 1.d0/AvN*1.d-3;
JMolToKcalMol = 0.0002390057361376673;
EhToKcalMol   = 627.5096080305927;


% %%% System Properties
% % O2+O 
% mO16      = 29148.94559d0; %15.9994d-3;
% mO18      = mO16 * 18.0 / 16.0;
% mO18O18   = mO18 + mO18;
% mu        = (mO16 * mO18O18) / (mO16 + mO18O18);
% 
% vO16_Lab    = 8026.0;
% vO18O18_Lab = 587.0;
% 
% Masses       = [mO18, mO18, mO16];



% %% CO+O
RunFldr  = '/home/venturi/WORKSPACE/CoarseAIR/CO2_ALL_SCATTERING/'
mO16     = 29148.94559d0;%15.9994d-3;
mO18     = mO16 * 18.0 / 16.0;
mC       = 21868.661757d0;%12.011e-3;
mCO18    = mC + mO18;
mu       = (mO16 * mCO18) / (mO16 + mCO18);

vO16_Lab    = 8100.0;
vCO18_Lab   = 800;

Masses       = [mC, mO16, mO16];



%%% Reading Quantum Mechanics
opts = delimitedTextImportOptions("NumVariables", 10);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["iTraj", "iPES", "Var3", "b_i", "Var5", "Var6", "Var7", "j_f", "v_f", "arr_f"];
opts.SelectedVariableNames = ["iTraj", "iPES", "b_i", "j_f", "v_f", "arr_f"];
opts.VariableTypes = ["double", "double", "string", "double", "string", "string", "string", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var3", "Var5", "Var6", "Var7"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var5", "Var6", "Var7"], "EmptyFieldRule", "auto");
tbl = readtable(strcat(RunFldr, "/Test/EMu/Bins_2_0/trajectories.csv"), opts);
Idx   = tbl.iTraj;
iPES  = tbl.iPES;
b_i   = tbl.b_i;
j_f   = tbl.j_f;
v_f   = tbl.v_f;
arr_f = tbl.arr_f;
clear opts tbl

NTrajsOI = min(length(arr_f), NTrajs);

Mapping = zeros(max(Idx),4);
iInel   = 0;
iExch1  = 0;
iExch2  = 0;
for iTraj = 1:NTrajsOI%length(Idx) 
    Mapping(Idx(iTraj),2) = iTraj; 
    Trajs.b(iTraj)        = b_i(iTraj);
    Trajs.iPES(iTraj)     = iPES(iTraj);
    Trajs.J(iTraj)        = j_f(iTraj);
    Trajs.v(iTraj)        = v_f(iTraj);
    
    if     (round(arr_f(iTraj)-0.5) == 16) || (round(arr_f(iTraj)-0.5) == 17)
        Mapping(Idx(iTraj),1) = 1; 
        iInel                 = iInel + 1; 
        Mapping(Idx(iTraj),3) = iInel; 
        
        Inelastic.b(iInel)    = b_i(iTraj);
        Inelastic.iPES(iInel) = iPES(iTraj);
        Inelastic.J(iInel)    = j_f(iTraj);
        Inelastic.v(iInel)    = v_f(iTraj);
    
    elseif (round(arr_f(iTraj)-0.5) == 32) || (round(arr_f(iTraj)-0.5) == 33)
        Mapping(Idx(iTraj),1) = 2;
        iExch1                = iExch1 + 1;
        Mapping(Idx(iTraj),3) = iExch1; 
        Mapping(Idx(iTraj),4) = iExch1+iExch2; 

        Exch1.b(iExch1)       = b_i(iTraj);
        Exch1.iPES(iExch1)    = iPES(iTraj);
        Exch1.J(iExch1)       = j_f(iTraj);
        Exch1.v(iExch1)       = v_f(iTraj);
        
        Exch.b(iExch1+iExch2)       = b_i(iTraj);
        Exch.iPES(iExch1+iExch2)    = iPES(iTraj);
        Exch.J(iExch1+iExch2)       = j_f(iTraj);
        Exch.v(iExch1+iExch2)       = v_f(iTraj);
    
    elseif (round(arr_f(iTraj)-0.5) == 48) || (round(arr_f(iTraj)-0.5) == 49)
        Mapping(Idx(iTraj),1) = 3;
        iExch2                = iExch2 + 1;
        Mapping(Idx(iTraj),3) = iExch2; 
        Mapping(Idx(iTraj),4) = iExch1+iExch2; 

        Exch2.b(iExch2)       = b_i(iTraj);        
        Exch2.iPES(iExch2)    = iPES(iTraj);
        Exch2.J(iExch2)       = j_f(iTraj);
        Exch2.v(iExch2)       = v_f(iTraj);

        Exch.b(iExch1+iExch2)       = b_i(iTraj);
        Exch.iPES(iExch1+iExch2)    = iPES(iTraj);
        Exch.J(iExch1+iExch2)       = j_f(iTraj);
        Exch.v(iExch1+iExch2)       = v_f(iTraj);
        
    end
    
end
NTrajs = max(Idx);
clear Idx

fprintf('Found %i    Inelastic Trajectories \n', iInel )
fprintf('Found %i 1st Exchange Trajectories \n', iExch1 )
fprintf('Found %i 2nd Exchange Trajectories \n', iExch2 )



%%% Reading Classical Mechanics
opts = delimitedTextImportOptions("NumVariables", 28);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Trajindex", "t_fin", "H_ini", "PaQ_ini1", "PaQ_ini2", "PaQ_ini3", "PaQ_ini4", "PaQ_ini5", "PaQ_ini6", "PaQ_ini7", "PaQ_ini8", "PaQ_ini9", "PaQ_ini10", "PaQ_ini11", "PaQ_ini12", "H_fin", "PaQ_fin1", "PaQ_fin2", "PaQ_fin3", "PaQ_fin4", "PaQ_fin5", "PaQ_fin6", "PaQ_fin7", "PaQ_fin8", "PaQ_fin9", "PaQ_fin10", "PaQ_fin11", "PaQ_fin12"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
tbl = readtable(strcat(RunFldr, "/Test/EMu/Bins_2_0/PaQSol.csv"), opts);
Idx       = tbl.Trajindex;
t_fin     = tbl.t_fin;
H_ini     = tbl.H_ini;
PaQ_ini1  = tbl.PaQ_ini1;
PaQ_ini2  = tbl.PaQ_ini2;
PaQ_ini3  = tbl.PaQ_ini3;
PaQ_ini4  = tbl.PaQ_ini4;
PaQ_ini5  = tbl.PaQ_ini5;
PaQ_ini6  = tbl.PaQ_ini6;
PaQ_ini7  = tbl.PaQ_ini7;
PaQ_ini8  = tbl.PaQ_ini8;
PaQ_ini9  = tbl.PaQ_ini9;
PaQ_ini10 = tbl.PaQ_ini10;
PaQ_ini11 = tbl.PaQ_ini11;
PaQ_ini12 = tbl.PaQ_ini12;
H_fin     = tbl.H_fin;
PaQ_fin1  = tbl.PaQ_fin1;
PaQ_fin2  = tbl.PaQ_fin2;
PaQ_fin3  = tbl.PaQ_fin3;
PaQ_fin4  = tbl.PaQ_fin4;
PaQ_fin5  = tbl.PaQ_fin5;
PaQ_fin6  = tbl.PaQ_fin6;
PaQ_fin7  = tbl.PaQ_fin7;
PaQ_fin8  = tbl.PaQ_fin8;
PaQ_fin9  = tbl.PaQ_fin9;
PaQ_fin10 = tbl.PaQ_fin10;
PaQ_fin11 = tbl.PaQ_fin11;
PaQ_fin12 = tbl.PaQ_fin12;
clear opts tbl
for iTraj = 1:NTrajsOI%length(Idx)
    if ( mod(iTraj,floor(NTrajsOI/20)) == 0 )
        fprintf('%f%% of Trajectories Postprocessed \n', iTraj/NTrajsOI*100 )
    end
    
    jTraj = Mapping(Idx(iTraj),2);
    if (jTraj ~= 0)

        Trajs.tf(jTraj)      = t_fin(iTraj);
        Trajs.Hi(jTraj)      = H_ini(iTraj);
        Trajs.Hf(jTraj)      = H_fin(iTraj);
        Trajs.PaQi(jTraj,1)  = PaQ_ini1(iTraj);
        Trajs.PaQi(jTraj,2)  = PaQ_ini2(iTraj);
        Trajs.PaQi(jTraj,3)  = PaQ_ini3(iTraj);
        Trajs.PaQi(jTraj,4)  = PaQ_ini4(iTraj);
        Trajs.PaQi(jTraj,5)  = PaQ_ini5(iTraj);
        Trajs.PaQi(jTraj,6)  = PaQ_ini6(iTraj);
        Trajs.PaQi(jTraj,7)  = PaQ_ini7(iTraj);
        Trajs.PaQi(jTraj,8)  = PaQ_ini8(iTraj);
        Trajs.PaQi(jTraj,9)  = PaQ_ini9(iTraj);
        Trajs.PaQi(jTraj,10) = PaQ_ini10(iTraj);
        Trajs.PaQi(jTraj,11) = PaQ_ini11(iTraj);
        Trajs.PaQi(jTraj,12) = PaQ_ini12(iTraj);
        Trajs.PaQf(jTraj,1)  = PaQ_fin1(iTraj);
        Trajs.PaQf(jTraj,2)  = PaQ_fin2(iTraj);
        Trajs.PaQf(jTraj,3)  = PaQ_fin3(iTraj);
        Trajs.PaQf(jTraj,4)  = PaQ_fin4(iTraj);
        Trajs.PaQf(jTraj,5)  = PaQ_fin5(iTraj);
        Trajs.PaQf(jTraj,6)  = PaQ_fin6(iTraj);
        Trajs.PaQf(jTraj,7)  = PaQ_fin7(iTraj);
        Trajs.PaQf(jTraj,8)  = PaQ_fin8(iTraj);
        Trajs.PaQf(jTraj,9)  = PaQ_fin9(iTraj);
        Trajs.PaQf(jTraj,10) = PaQ_fin10(iTraj);
        Trajs.PaQf(jTraj,11) = PaQ_fin11(iTraj);
        Trajs.PaQf(jTraj,12) = PaQ_fin12(iTraj);

        if     Mapping(Idx(iTraj),1) == 1
            iInel = Mapping(Idx(iTraj),3);

            Inelastic.tf(iInel)     = t_fin(iTraj);
            Inelastic.Hi(iInel)     = H_ini(iTraj);
            Inelastic.Hf(iInel)     = H_fin(iTraj);
            Inelastic.PaQi(iInel,:) = Trajs.PaQi(jTraj,:);
            Inelastic.PaQf(iInel,:) = Trajs.PaQf(jTraj,:);

            [X, Vi]                    = Transform_PaQ_To_XVCM(Inelastic.PaQi(iInel,:), Masses);
            Trajs.EColli_CM(jTraj)     = Compute_CollisionEnergy_CM(Vi, Masses, 1);
            Inelastic.EColli_CM(iInel) = Trajs.EColli_CM(jTraj);
            
            [X, Vf]                    = Transform_PaQ_To_XVCM(Inelastic.PaQf(iInel,:), Masses);
            Inelastic.ECollf_CM(iInel) = Compute_CollisionEnergy_CM(Vf, Masses, 1);
                        
            [Inelastic.Theta_CM(iInel), Inelastic.TOF(iInel)] = Compute_ScatteringAngle_CM(Vi, Vf, Masses, 1, Dist);
            
        elseif Mapping(Idx(iTraj),1) == 2
            iExch1 = Mapping(Idx(iTraj),3);
            iExch  = Mapping(Idx(iTraj),4);

            Exch1.tf(iExch1)     = t_fin(iTraj);
            Exch1.Hi(iExch1)     = H_ini(iTraj);
            Exch1.Hf(iExch1)     = H_fin(iTraj);
            Exch1.PaQi(iExch1,:) = Trajs.PaQi(jTraj,:);
            Exch1.PaQf(iExch1,:) = Trajs.PaQf(jTraj,:);
            
            [X, Vi]                 = Transform_PaQ_To_XVCM(Exch1.PaQi(iExch1,:), Masses);
            Trajs.EColli_CM(jTraj)  = Compute_CollisionEnergy_CM(Vi, Masses, 1);
            Exch1.EColli_CM(iExch1) = Trajs.EColli_CM(jTraj);
            
            [X, Vf]                 = Transform_PaQ_To_XVCM(Exch1.PaQf(iExch1,:), Masses);
            Exch1.ECollf_CM(iExch1) = Compute_CollisionEnergy_CM(Vf, Masses, 2);
                        
            [Exch1.Theta_CM(iExch1), Exch1.TOF(iExch1)] = Compute_ScatteringAngle_CM(Vi, Vf, Masses, 2, Dist);
            
            Exch.tf(iExch)        = Exch1.tf(iExch1);
            Exch.Hi(iExch)        = Exch1.Hi(iExch1);
            Exch.Hf(iExch)        = Exch1.Hf(iExch1);
            Exch.PaQi(iExch,:)    = Exch1.PaQi(iExch1,:);
            Exch.PaQf(iExch,:)    = Exch1.PaQf(iExch1,:);
            Exch.EColli_CM(iExch) = Exch1.EColli_CM(iExch1);
            Exch.ECollf_CM(iExch) = Exch1.ECollf_CM(iExch1);
            Exch.Theta_CM(iExch)  = Exch1.Theta_CM(iExch1);

        elseif Mapping(Idx(iTraj),1) == 3
            iExch2 = Mapping(Idx(iTraj),3);
            iExch  = Mapping(Idx(iTraj),4);

            Exch2.tf(iExch2)     = t_fin(iTraj);
            Exch2.Hi(iExch2)     = H_ini(iTraj);
            Exch2.Hf(iExch2)     = H_fin(iTraj);
            Exch2.PaQi(iExch2,:) = Trajs.PaQi(jTraj,:);
            Exch2.PaQf(iExch2,:) = Trajs.PaQf(jTraj,:);
            
            [X, Vi]                 = Transform_PaQ_To_XVCM(Exch2.PaQi(iExch2,:), Masses);
            Trajs.EColli_CM(jTraj)  = Compute_CollisionEnergy_CM(Vi, Masses, 1);
            Exch2.EColli_CM(iExch2) = Trajs.EColli_CM(jTraj);
            
            [X, Vf]                 = Transform_PaQ_To_XVCM(Exch2.PaQf(iExch2,:), Masses);
            Exch2.ECollf_CM(iExch2) = Compute_CollisionEnergy_CM(Vf, Masses, 3);
            
            [Exch2.Theta_CM(iExch2), Exch2.TOF(iExch2)] = Compute_ScatteringAngle_CM(Vi, Vf, Masses, 3, Dist);
            
            
            Exch.tf(iExch)        = Exch2.tf(iExch2);
            Exch.Hi(iExch)        = Exch2.Hi(iExch2);
            Exch.Hf(iExch)        = Exch2.Hf(iExch2);
            Exch.PaQi(iExch,:)    = Exch2.PaQi(iExch2,:);
            Exch.PaQf(iExch,:)    = Exch2.PaQf(iExch2,:);
            Exch.EColli_CM(iExch) = Exch2.EColli_CM(iExch2);
            Exch.ECollf_CM(iExch) = Exch2.ECollf_CM(iExch2);
            Exch.Theta_CM(iExch)  = Exch2.Theta_CM(iExch2);
            
        end
        
    end
    
end



figure(1)
histogram(Trajs.EColli_CM, 100)
hold on

figure(2)
histogram(Inelastic.ECollf_CM, 100)
hold on

figure(3)
histogram(Exch.ECollf_CM, 100)
hold on
%histogram(Exch2.ECollf_CM, 100)


figure(4)
histogram(Inelastic.Theta_CM, 100)
hold on

figure(5)
histogram(Exch.Theta_CM, 100)
hold on
%histogram(Exch2.Theta_CM, 100)




figure(6)
histogram(Inelastic.b, 100)
hold on


figure(7)
histogram(Exch.b, 100)
hold on
%histogram(Exch2.b, 100)



bVec      = linspace(0.0, max(Trajs.b), 50);
bVec_Inel = histcounts(Inelastic.b, bVec);
bVec_Exch = histcounts(Exch.b, bVec);
bVec_Tot  = bVec_Inel + bVec_Exch;
figure(8)
semilogy(bVec(1:end-1), bVec_Exch./bVec_Tot)
hold on
%histogram(Exch2.b, 100)


figure(10)
scatter(Exch.b, Exch.Theta_CM)

figure(11)
scatter(Exch.EColli_CM, Exch.Theta_CM)