close all
clear all
clc


RunFldr  = "/home/venturi/Desktop/";
TT       = 15000;
Velocity = 0.1247967544D-10;


opts = delimitedTextImportOptions("NumVariables", 14);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["iTraj", "iPES", "bmax", "b_i", "j1_i", "v1_i", "j2_i", "v2_i", "arr_i", "j1_f", "v1_f", "j2_f", "v2_f", "arr_f"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
tbl = readtable(strcat(RunFldr,"/trajectories.csv"), opts);
Idx   = tbl.iTraj;
iPES  = tbl.iPES;
bmax  = tbl.bmax;
b_i   = tbl.b_i;
j1_i  = tbl.j1_i;
v1_i  = tbl.v1_i;
j2_i  = tbl.j2_i;
v2_i  = tbl.v2_i;
arr_i = tbl.arr_i;
j1_f  = tbl.j1_f;
v1_f  = tbl.v1_f;
j2_f  = tbl.j2_f;
v2_f  = tbl.v2_f;
arr_f = tbl.arr_f;
clear opts tbl
NTrajs  = length(Idx);
bUnique = unique(bmax)';
Nbb     = length(bUnique);

for ib=1:Nbb
   TrajVec(ib).bMax      = bUnique(ib);
   TrajVec(ib).ACircle   = pi * bUnique(ib)^2;
   TrajVec(ib).TotNb     = 0;
   TrajVecTemp(ib).b_i   = zeros(NTrajs,1);
   TrajVecTemp(ib).j1_f  = zeros(NTrajs,1);
   TrajVecTemp(ib).v1_f  = zeros(NTrajs,1);
   TrajVecTemp(ib).j2_f  = zeros(NTrajs,1);
   TrajVecTemp(ib).v2_f  = zeros(NTrajs,1);
   TrajVecTemp(ib).arr_f = zeros(NTrajs,1);
end
TrajVec(1).ARing = TrajVec(1).ACircle;
for ib=2:Nbb
   TrajVec(ib).ARing = TrajVec(ib).ACircle - TrajVec(ib-1).ACircle;
end

for iTraj=1:NTrajs
   ib                = sum(b_i(iTraj) > bUnique) + 1;
   if (not(isnan(arr_f(iTraj))))
       TrajVec(ib).TotNb = TrajVec(ib).TotNb + 1;
       TrajVecTemp(ib).b_i(TrajVec(ib).TotNb)   = b_i(iTraj);
       TrajVecTemp(ib).j1_f(TrajVec(ib).TotNb)  = j1_f(iTraj);
       TrajVecTemp(ib).v1_f(TrajVec(ib).TotNb)  = v1_f(iTraj);
       TrajVecTemp(ib).j2_f(TrajVec(ib).TotNb)  = j2_f(iTraj);
       TrajVecTemp(ib).v2_f(TrajVec(ib).TotNb)  = v2_f(iTraj);
       TrajVecTemp(ib).arr_f(TrajVec(ib).TotNb) = arr_f(iTraj);
   end
end
for ib=1:Nbb
   TrajVec(ib).b_i     = TrajVecTemp(ib).b_i(1:TrajVec(ib).TotNb);
   TrajVec(ib).j1_f    = TrajVecTemp(ib).j1_f(1:TrajVec(ib).TotNb);
   TrajVec(ib).v1_f    = TrajVecTemp(ib).v1_f(1:TrajVec(ib).TotNb);
   TrajVec(ib).j2_f    = TrajVecTemp(ib).j2_f(1:TrajVec(ib).TotNb);
   TrajVec(ib).v2_f    = TrajVecTemp(ib).v2_f(1:TrajVec(ib).TotNb);
   TrajVec(ib).arr_f   = TrajVecTemp(ib).arr_f(1:TrajVec(ib).TotNb);
end
clear TrajVecTemp

for iP=1:3
    Pair(iP).Type_vs_b     = [];
    Pair(iP).CrossSec_vs_b = [];
end
for ib=1:Nbb
   for iP=1:3
       TrajVec(ib).Pair(iP).TotNb = 0;
       TrajVec(ib).Pair(iP).TypeTemp = zeros(TrajVec(ib).TotNb,1);
   end
   TrajVec(ib).iP  = floor(TrajVec(ib).arr_f/16);
   
   for iTraj=1:TrajVec(ib).TotNb
       iP                               = TrajVec(ib).iP(iTraj);
       arr_f                            = TrajVec(ib).arr_f(iTraj);
       TypeA                            = ceil( ceil(mod(mod(arr_f,16),4))   / 2);
       TypeB                            = ceil( (floor(mod(arr_f,16)/4) + 1) / 2);
       if TypeA < 1 || TypeA > 2 || TypeB < 1 || TypeB > 2
           pause
       end
       TrajVec(ib).Pair(iP).TotNb           = TrajVec(ib).Pair(iP).TotNb + 1;
       jTraj                                = TrajVec(ib).Pair(iP).TotNb;
       TrajVec(ib).Pair(iP).TypeTemp(jTraj) = (TypeA-1).*2 + TypeB;
   end 
   
   for iP=1:3
        TrajVec(ib).Pair(iP).Type     = TrajVec(ib).Pair(iP).TypeTemp(1:TrajVec(ib).Pair(iP).TotNb);
        
        TrajVec(ib).Pair(iP).Count    = histcounts(TrajVec(ib).Pair(iP).Type, [0.5,1.5,2.5,3.5,4.5]) ./ max(1.e-10,TrajVec(ib).TotNb);
        Pair(iP).Type_vs_b            = [Pair(iP).Type_vs_b;     TrajVec(ib).Pair(iP).Count];

        TrajVec(ib).Pair(iP).CrossSec = TrajVec(ib).Pair(iP).Count .* TrajVec(ib).ARing;
        Pair(iP).CrossSec_vs_b        = [Pair(iP).CrossSec_vs_b; TrajVec(ib).Pair(iP).CrossSec];
   end
%    figure(1)
%    histogram(TrajVec(ib).Pair(1).Type, 'Normalization','pdf')
%    hold on
%    histogram(TrajVec(ib).Pair(2).Type, 'Normalization','pdf')
%    histogram(TrajVec(ib).Pair(3).Type, 'Normalization','pdf')
%    hold off
%    pause
   
end
for iP=1:3
    Pair(iP).CrossSec = sum(Pair(iP).CrossSec_vs_b,1);
    Pair(iP).Rates    = Pair(iP).CrossSec .* Velocity;
end

iP=1;
figure((iP-1)*10+1)
bar(Pair(iP).Type_vs_b)
legend('Inelastic','Simple','Simple','Double')

figure((iP-1)*10+2)
bar(Pair(iP).CrossSec_vs_b)
legend('Inelastic','Simple','Simple','Double')

iP=2;
figure((iP-1)*10+1)
bar(Pair(iP).Type_vs_b)
legend('Inelastic','Swap 1, Type 1','Swap 2, Type 1','Double')

figure((iP-1)*10+2)
bar(Pair(iP).CrossSec_vs_b)
legend('Inelastic','Swap 1, Type 1','Swap 2, Type 1','Double')

iP=3;
figure((iP-1)*10+1)
bar(Pair(iP).Type_vs_b)
legend('Inelastic','Swap 1, Type 2','Swap 2, Type 2','Double')

figure((iP-1)*10+2)
bar(Pair(iP).CrossSec_vs_b)
legend('Inelastic','Swap 1, Type 2','Swap 2, Type 2','Double')


RatesSimple = Pair(1).Rates(2)
RatesSwap   = Pair(2).Rates(2) + Pair(3).Rates(2)
RatesDouble = (Pair(1).Rates(4) + Pair(2).Rates(4) + Pair(3).Rates(4)) / 2
RatesTot    = RatesSimple + RatesSwap + RatesDouble

PercSimple  = RatesSimple / RatesTot * 100.0
PercSwap    = RatesSwap   / RatesTot * 100.0
PercDouble  = RatesDouble / RatesTot * 100.0







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison with Ross
%%

TVec          = [4000, 5000, 6500, 8000, 10000, 13000];
TSpace        = linspace(4000,16000,100);
KDVec         = [2.79e-15, 4.78e-14, 6.40e-13, 3.163e-12, 1.2147e-11, 3.948e-11];
PercSimpleVec = [86.4, 84.0, 80.38, 76.59, 72.59, 67.63];
PercSwapVec   = [13.6, 16.0, 19.60, 23.32, 27.00, 31.07];
PercDoubleVec = [0.0, 0.0, 0.021, 0.095, 0.409, 1.301];

[xData, yData] = prepareCurveData( 10000./TVec, log10(KDVec) );
ft = fittype( 'pchipinterp' );
[fitresultKD, gof] = fit( xData, yData, ft );

[xData, yData] = prepareCurveData( TVec, PercSimpleVec );
ft = fittype( 'pchipinterp' );
[fitresultSimple, gof] = fit( xData, yData, ft );

[xData, yData] = prepareCurveData( TVec, PercSwapVec );
ft = fittype( 'pchipinterp' );
[fitresultSwap, gof] = fit( xData, yData, ft );

[xData, yData] = prepareCurveData( TVec, PercDoubleVec );
ft = fittype( 'pchipinterp' );
[fitresultDouble, gof] = fit( xData, yData, ft );


figure(101)
semilogy(10000./TVec, KDVec, '.');
hold on
semilogy(10000./TSpace, 10.0.^fitresultKD(10000./TSpace), '-');
semilogy(10000./TT,   RatesTot, 'o-');


figure(102)
plot(TVec, PercSimpleVec,'.k');
hold on
plot(TVec, PercSwapVec,'.r');
plot(TVec, PercDoubleVec,'.b');

plot(TSpace, fitresultSimple(TSpace),'-k');
hold on
plot(TSpace, fitresultSwap(TSpace),'-r');
plot(TSpace, fitresultDouble(TSpace),'-b');

plot(TT, PercSimple,'o-k');
plot(TT, PercSwap,'o-r');
plot(TT, PercDouble,'o-b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%