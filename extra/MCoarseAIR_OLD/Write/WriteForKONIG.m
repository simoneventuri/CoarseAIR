clc
%close all


iBinnedMol = 1;
iFigure    = 2001;
iSteps     = 1;


clear Rates
Rates(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol)) = ( RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),1) + RatesMatrix(1:NLevels(iBinnedMol),1:NLevels(iBinnedMol),2) )  .* 1.d-6;
DissRates                                          = ProcessesRates(1:NLevels(iBinnedMol),1,1)' .* 1.d-6;


LevToBinFinal = LevToBin(:,iBinnedMol,iSteps);
NbBins        = max(LevToBinFinal);

ExpVec0(1:NLevels(iBinnedMol),iBinnedMol) = exp( - LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (300.d0 .* UKb) );
Q0(1:NLevels(iBinnedMol),iBinnedMol)      = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* ExpVec0(1:NLevels(iBinnedMol),iBinnedMol);

ExpVec(1:NLevels(iBinnedMol),iBinnedMol)  = exp( - LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol) .* Ue ./ (T0_Vec(1) .* UKb) );
Q(1:NLevels(iBinnedMol),iBinnedMol)       = Levelg(1:NLevels(iBinnedMol),iBinnedMol) .* ExpVec(1:NLevels(iBinnedMol),iBinnedMol);
QMat                                      = kron(Q,1.d0./Q');  
QEn(1:NLevels(iBinnedMol),iBinnedMol)     = Q(1:NLevels(iBinnedMol),iBinnedMol) .* LevelEeV0(1:NLevels(iBinnedMol),iBinnedMol);

QBin  = zeros(NbBins,1);
QBin0 = zeros(NbBins,1);
EeVBin = zeros(NbBins,1);
for iLevels = 1:NLevels(iBinnedMol)
  QBin(LevToBinFinal(iLevels))   = QBin(LevToBinFinal(iLevels))   + Q(iLevels);
  QBin0(LevToBinFinal(iLevels))  = QBin0(LevToBinFinal(iLevels))  + Q0(iLevels);
  EeVBin(LevToBinFinal(iLevels)) = EeVBin(LevToBinFinal(iLevels)) + QEn(iLevels);
end

BinsRates     = zeros(NbBins,NbBins);
BinsDissRates = zeros(NbBins,1);
for i = 1:size(Rates,1)
  for j = 1:size(Rates,2)
    if LevToBinFinal(i) ~= LevToBinFinal(j)
      if i > j
        BinsRates(LevToBinFinal(i),LevToBinFinal(j)) = BinsRates(LevToBinFinal(i),LevToBinFinal(j)) + (Q(i,iBinnedMol) ./ QBin(LevToBinFinal(i))) .* Rates(i,j);
      else
        BinsRates(LevToBinFinal(i),LevToBinFinal(j)) = BinsRates(LevToBinFinal(i),LevToBinFinal(j)) + (Q(i,iBinnedMol) ./ QBin(LevToBinFinal(i))) .* Rates(j,i) .* QMat(j,i);
      end
    end
  end
  BinsDissRates(LevToBinFinal(i)) = BinsDissRates(LevToBinFinal(i)) + (Q(i,iBinnedMol) ./ QBin(LevToBinFinal(i))) .* DissRates(i);
end



FileFormatName = strcat('/Users/sventuri/WORKSPACE/neqplasma_QCT/database_CO2Binned/thermo/CO_Format')
FileFormatID   = fopen(FileFormatName);

FileNewName    = strcat('/Users/sventuri/WORKSPACE/neqplasma_QCT/database_CO2Binned/thermo/CO')
fileNewID      = fopen(FileNewName,'w');

tline  = fgetl(FileFormatID);
ttline = strcat(tline,'\n');
while ischar(tline)
  fprintf(fileNewID,ttline);
  tline = fgetl(FileFormatID);
  ttline = strcat(tline,'\n');
end
fclose(FileFormatID);

fprintf(fileNewID,'NB_ENERGY_LEVELS = %g\n', NbBins);
for iBin = 1:NbBins
  fprintf(fileNewID,'%14.8E  %14.8E\n', QBin(iBin), 0.d0);%EeVBin(iBin))
end
fclose(fileNewID);  



FileNewName    = strcat('/Users/sventuri/WORKSPACE/neqplasma_QCT/database_CO2Binned/kinetics/CO2')
fileNewID      = fopen(FileNewName,'w');

fprintf(fileNewID,'Units=cm^3/s\n');
for iBin = 1:NbBins
  for jBin = 1:iBin-1
    if BinsRates(iBin,jBin) > 0.d0
      fprintf(fileNewID,'CO(%g)+O=CO(%g)+O:%+11.4E,%+11.4E,%+11.4E,6\n', iBin, jBin, BinsRates(iBin,jBin) * 1.d6, 0.d0, 0.d0);
    end
  end
  %if BinsDissRates(iBin) > 0.d0
  %  fprintf(fileNewID,'CO(%g)+O=C+O+O:%+11.4E,%+11.4E,%+11.4E,2\n', iBin, BinsDissRates(iBin) * 1.d6, 0.d0, 0.d0);
  %end
end
fclose(fileNewID);  



FileNewName    = strcat('/Users/sventuri/WORKSPACE/neqplasma_QCT/konig/runs/CO2/input/InitialMolFracs.dat')
fileNewID      = fopen(FileNewName,'w');

fprintf(fileNewID,'Nb of Components:\n');
fprintf(fileNewID,'  3\n');
fprintf(fileNewID,'Nb of Levels per Component:\n');
fprintf(fileNewID,'  1\n');
fprintf(fileNewID,'  1\n');
fprintf(fileNewID,'  %g\n', NbBins);
fprintf(fileNewID,'Nb of Levels per Component:\n');
fprintf(fileNewID,'  %14.8E\n', 1.d0);
fprintf(fileNewID,'  %14.8E\n', 1.d0);
for iBin = 1:NbBins
  fprintf(fileNewID,'  %14.8E\n', QBin0(iBin) ./ sum(QBin0) );
end
fclose(fileNewID);  