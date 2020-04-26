function [tStart, iQSS, tEnd] =  FindQSS(iT, t, ProcessesRatesOverall)

  global T0_Vec
  global DissFlg InelFlg ExchFlg

%   if T0_Vec < 5000
%     KQSS_Eps = 1.e-9;
%     EpsT     = 1.e-5;
%     EpsTT    = 6.e-4;
%     EpsTTT   = 2.e-3;
%   elseif T0_Vec < 12000
%     KQSS_Eps = 1.e-9;
%     EpsT     = 3.e-5;
%     EpsTT    = 6.e-4;
%     EpsTTT   = 2.e-3;
%   else
%     KQSS_Eps = 1.e-9;
%     EpsT     = 5.e-5;
%     EpsTT    = 6.e-4;
%     EpsTTT   = 2.e-3;
%   end
%   
  
  KQSS_Eps = 1.e-9;
  EpsT     = 5.e-4;%5.e-8;
  EpsTT    = 1.e-3;
  InpPerc  = 0.1;
  
  yy   = ProcessesRatesOverall(:,1);

  
  iQSS = 1;
  while (yy(iQSS) <= yy(1) + 0.1*(yy(end) - yy(1)))
      iQSS = iQSS + 1;
  end
  while ( abs(log10(yy(iQSS)) - log10(yy(iQSS+1))) / abs(log10(yy(iQSS))) > EpsT ) && ( iQSS < length(yy)-3)
      iQSS = iQSS + 1;
  end
  iQSS_Start = iQSS - 1;
  while ( abs(log10(yy(iQSS)) - log10(yy(iQSS+1))) / abs(log10(yy(iQSS))) <= EpsTT ) && ( iQSS < length(yy)-3)
      iQSS = iQSS + 1;
  end
  iQSS_End = iQSS - 1;
  iQSS     = floor( (iQSS_Start + iQSS_End) / 2.0 );
  
  clear fitresult yyy
  SStart = 1;
  while (yy(SStart) < yy(iQSS)/6)
    SStart = SStart + 1;
  end
  EEnd = size(yy,1);%iQSS;
%   while (abs(yy(EEnd) - yy(iQSS))/yy(iQSS) < (1.0+InpPerc+0.05))
%     EEnd = EEnd + 1;
%   end
  yyy(SStart:EEnd)  = yy(SStart:EEnd) - yy(iQSS)*(1.0-InpPerc);
  [xData, yData] = prepareCurveData( t(SStart:EEnd), yyy(SStart:EEnd) );
  ft = 'splineinterp';
  [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' )
  tStart = fzero(fitresult, t(iQSS)/3.0)
  
  clear fitresult yyy xData yData
  yyy(SStart:EEnd)  = yy(SStart:EEnd) - yy(iQSS)*(1.0+InpPerc);
  [xData, yData] = prepareCurveData( t(SStart:EEnd), yyy(SStart:EEnd) );
  ft = 'splineinterp';
  [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' )
  tEnd   = fzero(fitresult, t(iQSS)*5)

% 
%   jQSS = iQSS_Start;
%   while abs(log10(yy(jQSS)) - log10(yy(iQSS))) > 0.05
%       jQSS = jQSS + 1;
%   end
%   iQSS_Start = jQSS;
%   jQSS       = iQSS;
%   while abs( log10(yy(jQSS)) - log10(yy(iQSS)) ) < 0.05
%       jQSS = jQSS + 1;
%   end
%   iQSS_End   = jQSS;
%   
  
  FileName = strcat('./KAveraged_',num2str(DissFlg),'_',num2str(InelFlg),'_',num2str(ExchFlg),'_',num2str(ExchFlg),'.csv');
  if exist(FileName, 'file')
    fileID1  = fopen(FileName,'a');
  else
    fileID1  = fopen(FileName,'w');
    fprintf(fileID1,'# T [K], KDiss_Eq, KExch_Eq, KDiss_QSS, KExch_QSS\n');
  end
  fprintf(fileID1,'%e,%e,%e,%e,%e\n', T0_Vec(1), ProcessesRatesOverall(end,1,1), ProcessesRatesOverall(end,3,1)+ProcessesRatesOverall(end,4,1), ProcessesRatesOverall(iQSS,1,1), ProcessesRatesOverall(iQSS,3,1)+ProcessesRatesOverall(iQSS,4,1) )
  fclose(fileID1);
  
  
  figure; 
  semilogx(t,log10(yy))
  hold on
  semilogx(tStart,log10(yy(iQSS)*(1.0-InpPerc)),'o')
  semilogx(tEnd,log10(yy(iQSS)*(1.0+InpPerc)),'o')
  semilogx(t(iQSS),log10(yy(iQSS)),'o')
  
end