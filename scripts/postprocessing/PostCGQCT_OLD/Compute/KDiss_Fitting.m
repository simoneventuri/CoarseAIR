close all
clc

global T0_Vec System SystemPath 
global NBins BinnedMolName NBinnedMol KeV KinMthd EhToeV MoleculesName AtomMass Pair_to_Atoms ParaViewFlg


%% Excluding Levels
RejectedLevels = [];


%% Defining Input Variables 
clear x1 x2 y
x1     = log10(-DeltaEintDiss); %log10(abs(LevelEeV - VMax));
x1Name = 'Energy Distance from J Barrier [eV]';

x2     = rEquival;%(rOut + rIn)./2 ./ Tau; %rOut; %rEquival; %(rOut + rIn)./2; %./ Tau;
x2Name = 'Equivalent R [a_0]';
%x2 = log10(Tau);


%% Defining Output Variables
KDiss     = ProcessesRates(:,1);
LogKDiss  = log10(ProcessesRates(:,1));
y         = LogKDiss;
yName     = 'Log(KDiss_i)';
KDiss_Min = 1.d-14;

FileName = strcat('./',BinnedMolName,'_KDiss_Orig.csv');
fileID   = fopen(FileName,'w');
TempStr  =  strcat(BinnedMolName,'(%i)+',AtomsName(1),'=',AtomsName(1),'+',AtomsName(2),'+',AtomsName(3),':+%e,+0.0000E+00,+0.0000E+00,2\n');
for iLevels=1:size(Leveljqn,1)
  if (KDiss(iLevels) > KDiss_Min)
    fprintf(fileID, TempStr, iLevels, KDiss(iLevels));
  end
end
fclose(fileID);


%% Plotting Variables
figure(101)
hold on
scatter(x1, y, 'bo', 'filled')
xlabel(x1Name)
ylabel(yName)

figure(102)
hold on
scatter(x2, y, 'ro', 'filled')
xlabel(x2Name)
ylabel(yName)




%% Removing Levels
clear x1_New x2_New y_New
j=0;
for iLevels=1:length(LogKDiss)
   if (y(iLevels) >= log10(KDiss_Min)) && (sum(iLevels ~= RejectedLevels) == 0)
    j=j+1;
    x1_New(j) = x1(iLevels);
    x2_New(j) = x2(iLevels);
    y_New(j)  = y(iLevels);
  end
end


%% Fitting and Predicting
f     = fit( [x1_New', x2_New'], y_New', 'poly22' );
yPred = f.p00        + ...
        f.p10.*x1    + f.p01.*x2    + ...
        f.p20.*x1.^2 + f.p02.*x2.^2 + f.p11.*x1.*x2;%    + ...
        %f.p30.*x1.^3 + f.p03.*x2.^3 + f.p21.*x1.^2.*x2 + f.p12.*x1.*x2.^2 + ...
        %f.p40.*x1.^4 + f.p04.*x2.^4 + f.p31.*x1.^3.*x2 + f.p13.*x1.*x2.^3 + f.p22.*x1.^2.*x2.^2;% + ...
        %f.p50.*x1.^5 + f.p05.*x2.^5 + f.p41.*x1.^4.*x2 + f.p14.*x1.*x2.^4 + f.p32.*x1.^3.*x2.^2 + f.p23.*x1.^2.*x2.^3; % + ...
        %f.p60.*x1.^6 + f.p06.*x2.^6 + f.p51.*x1.^5.*x2 + f.p15.*x1.*x2.^5 + f.p42.*x1.^4.*x2.^2 + f.p24.*x1.^2.*x2.^4 + f.p33.*x1.^3.*x2.^3 );
yPredName = 'Predicted Log(KDiss_i)';
KDissPred = 10.d0.^yPred;

% Plotting Scatter Plot                    
figure(201)
scatter(y,yPred,'o', 'filled')
hold on
plot([log10(KDiss_Min),-8],[log10(KDiss_Min),-8],'-')
xlim([log10(KDiss_Min),-8])
ylim([log10(KDiss_Min),-8])
xlabel(yName)
ylabel(yPredName)

figure(202)
scatter(KDiss, KDissPred, 'o', 'filled')
hold on
plot([KDiss_Min,10^(-8)],1/2.*[KDiss_Min,10^(-8)],'r-')
plot([KDiss_Min,10^(-8)],2.*[KDiss_Min,10^(-8)],'r-')
plot([KDiss_Min,10^(-8)],[KDiss_Min,10^(-8)],'k-')
xlabel('KDiss_i')
ylabel('Predicted KDiss_i')


%% Writing Predicted Rates and Errors
FileName  = strcat('./',BinnedMolName,'_KDiss_Pred.csv');
fileID    = fopen(FileName,'w');
TempStr   =  strcat(BinnedMolName,'(%i)+',AtomsName(1),'=',AtomsName(1),'+',AtomsName(2),'+',AtomsName(3),':+%e,+0.0000E+00,+0.0000E+00,2\n');
for iLevels=1:size(Leveljqn,1)
  if (KDissPred(iLevels) > KDiss_Min)
    fprintf(fileID, TempStr, iLevels, KDissPred(iLevels));
  end
end
fclose(fileID);

KDissError = abs(KDissPred - KDiss)
FileName   = strcat('./',BinnedMolName,'_KDiss_Error.csv');
fileID     = fopen(FileName,'w');
fprintf(fileID,'#rIn,rOut,J,EeV,AbsError,Error\n');
for iLevels=1:size(Leveljqn,1)
  fprintf(fileID,'%e,%e,%e,%e,%e,%e\n', rIn(iLevels), rOut(iLevels), Leveljqn(iLevels), LevelEeV(iLevels), KDissError(iLevels)/KDiss(iLevels), KDissError(iLevels) );
end
fclose(fileID);


%% Plotting Rate Comparisons
figure(203)
scatter(LevelEeV, KDiss, 'ko', 'filled')
hold on
scatter(LevelEeV, KDissPred, 'bo', 'filled')
%scatter(LevelEeV, KDissError, 'ro', 'filled')
xlabel('Energy [eV]')
ylabel('KDiss')
set(gca,'yscale','log')




%% Selecting Set of Sampled Points
figure(301)
plot(x1,x2,'o')

[aSorted, aOrder] = sort(DeltaEintDiss);
x1Bound       = log10(abs([-1.e-200; -1.e-5; -1.e-3; -1.e-1; -0.5; -1.0; -2.2; -3.5; -5.0]));
NPoints       = length(x1Bound);
FinalId       = [];
for i=2:NPoints
  
  k=1;
  for j=1:length(x1)
    if (x1(j) >= x1Bound(i-1)) && (x1(j) < x1Bound(i))
      Tempx1(k) = x1(j);
      Tempx2(k) = x2(j);
      TempId(k) = j;
      k=k+1;
    end
  end
  
  NTemp             = length(Tempx2);
  NPoints2          = floor(NTemp/NPoints);
  [bSorted, bOrder] = sort(Tempx2);
  FinalId           = [FinalId, TempId(bOrder([1:NPoints2:NTemp]))];
  
  clear TempDeltaE TempDeltar TempId
end

x1_New2 = x1(FinalId);
x2_New2 = x2(FinalId);
y_New2  = y(FinalId);

hold on
plot(x1_New2,x2_New2,'ro')


%% Fitting and Predicting
g = fit( [x1_New2, x2_New2], y_New2, 'poly55' );

yPred2 = g.p00        + ...
         g.p10.*x1    + g.p01.*x2    + ...
         g.p20.*x1.^2 + g.p02.*x2.^2 + g.p11.*x1.*x2    + ...
         g.p30.*x1.^3 + g.p03.*x2.^3 + g.p21.*x1.^2.*x2 + g.p12.*x1.*x2.^2 + ...
         g.p40.*x1.^4 + g.p04.*x2.^4 + g.p31.*x1.^3.*x2 + g.p13.*x1.*x2.^3 + g.p22.*x1.^2.*x2.^2 + ...
         g.p50.*x1.^5 + g.p05.*x2.^5 + g.p41.*x1.^4.*x2 + g.p14.*x1.*x2.^4 + g.p32.*x1.^3.*x2.^2 + g.p23.*x1.^2.*x2.^3;% + ...
         %g.p60.*x1.^6 + g.p06.*x2.^6 + g.p51.*x1.^5.*x2 + g.p15.*x1.*x2.^5 + g.p42.*x1.^4.*x2.^2 + g.p24.*x1.^2.*x2.^4 + g.p33.*x1.^3.*x2.^3 );
KDissPred2 = 10.d0.^yPred2;

% Plotting Scatter Plot                    
figure(302)
scatter(y, yPred2, 'o', 'filled')
hold on
plot([log10(KDiss_Min),-8],[log10(KDiss_Min),-8],'-')
xlim([log10(KDiss_Min),-8])
ylim([log10(KDiss_Min),-8])
xlabel(yName)
ylabel(yPredName)

figure(303)
scatter(KDiss, KDissPred2, 'o', 'filled')
hold on
plot([KDiss_Min,10^(-8)],1/2.*[KDiss_Min,10^(-8)],'r-')
plot([KDiss_Min,10^(-8)],2.*[KDiss_Min,10^(-8)],'r-')
plot([KDiss_Min,10^(-8)],[KDiss_Min,10^(-8)],'k-')
xlabel('KDiss_i')
ylabel('Predicted KDiss_i')


%% Plotting Rate Comparisons
figure(304)
scatter(LevelEeV, KDiss, 'ko', 'filled')
hold on
scatter(LevelEeV, KDissPred2, 'ro', 'filled')
%scatter(LevelEeV, KDissError, 'ro', 'filled')
xlabel('Energy [eV]')
ylabel('KDiss')
set(gca,'yscale','log')