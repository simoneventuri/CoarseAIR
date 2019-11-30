close all
clear all
clc


WORKSPACE_PATH = '/Users/sventuri/WORKSPACE/';
PES_NAMES      = [string('11A1')];
NPESs          = size(PES_NAMES,2);
iFigure        = 1;

AngToBohr   = 1.d0 / 0.52917721092;
KcalMolToEV = 0.0433641153087705;
ConvConst   = 120.311d0;%120.1745883646d0;%120.31146333%120.1745883646%120.31146333%94.21146333;%94.0745884; 

% figure
iPES = 0;
for PES_NAME = PES_NAMES
  iPES = iPES + 1
  PES_NAME
  
  filename = strcat(WORKSPACE_PATH,'/CoarseAIR/coarseair/dtb/N2O2/PESs/AbInitioData/', PES_NAME(:), '.dat')
  delimiter = ' ';
  formatSpec = '%f%f%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
  fclose(fileID);
  Col1 = dataArray{:, 1};
  Col2 = dataArray{:, 2};
  Col3 = dataArray{:, 3};
  clearvars filename delimiter formatSpec fileID dataArray ans;
  for i=1:5:length(Col1)
    E_1(Col1(i)) = (Col2(i) - ConvConst).*KcalMolToEV;
    E_2(Col1(i)) = (Col3(i) - ConvConst).*KcalMolToEV;
  end
%   figure(iFigure)
%   scatter([1:size(E_1,2)],E_1)
%   iFigure = iFigure + 1;
%   pause

  iPoint=1;
  for i=2:5:length(Col1)
    X1(iPoint,:) = [Col1(i), Col2(i), Col3(i)];
    iPoint=iPoint+1;
  end

  iPoint=1;
  for i=3:5:length(Col1)
    X2(iPoint,:) = [Col1(i), Col2(i), Col3(i)];
    iPoint=iPoint+1;
  end

  iPoint=1;
  for i=4:5:length(Col1)
    X3(iPoint,:) = [Col1(i), Col2(i), Col3(i)];
    iPoint=iPoint+1;
  end
  
  iPoint=1;
  for i=5:5:length(Col1)
    X4(iPoint,:) = [Col1(i), Col2(i), Col3(i)];
    iPoint=iPoint+1;
  end
  
  R1(:) = sqrt( (X2(:,1)-X1(:,1)).^2 + (X2(:,2)-X1(:,2)).^2 + (X2(:,3)-X1(:,3)).^2 ) .* AngToBohr;
  R2(:) = sqrt( (X3(:,1)-X1(:,1)).^2 + (X3(:,2)-X1(:,2)).^2 + (X3(:,3)-X1(:,3)).^2 ) .* AngToBohr;
  R3(:) = sqrt( (X4(:,1)-X1(:,1)).^2 + (X4(:,2)-X1(:,2)).^2 + (X4(:,3)-X1(:,3)).^2 ) .* AngToBohr;
  R4(:) = sqrt( (X3(:,1)-X2(:,1)).^2 + (X3(:,2)-X2(:,2)).^2 + (X3(:,3)-X2(:,3)).^2 ) .* AngToBohr;
  R5(:) = sqrt( (X4(:,1)-X2(:,1)).^2 + (X4(:,2)-X2(:,2)).^2 + (X4(:,3)-X2(:,3)).^2 ) .* AngToBohr;
  R6(:) = sqrt( (X4(:,1)-X3(:,1)).^2 + (X4(:,2)-X3(:,2)).^2 + (X4(:,3)-X3(:,3)).^2 ) .* AngToBohr;
  
  FileName1 = strcat('./N2O_AbInitioPoints.csv');
  fileID1   = fopen(FileName1,'w');
  fprintf(fileID1,'R_N2,R_O2,R_O2,E1,E2\n');
  
  iDaje  = 0;
  MinVal = 15.0;
  for i=1:length(R1)
    if (R1(i) >= MinVal) && (R2(i) >= MinVal) && (R3(i) >= MinVal)
      iDaje=iDaje+1;
      fprintf(fileID1,'%e,%e,%e,%e,%e\n', R6(i), R5(i), R4(i), E_1(i), E_2(i));
    end
  end
  iDaje
  jDaje  = 0;
  for i=1:length(R1)
    if (R1(i) >= MinVal) && (R5(i) >= MinVal) && (R4(i) >= MinVal)
      jDaje=jDaje+1;
      fprintf(fileID1,'%e,%e,%e,%e,%e\n', R6(i), R2(i), R3(i), E_1(i), E_2(i));
    end
  end
  jDaje
  fclose(fileID1);
  pause
  
  
  R       = [R1', R2', R3', R4', R5', R6'];
  RSorted = sort(R,2);
  alpha(:,3) = acos( (RSorted(:,1).^2 + RSorted(:,2).^2 - RSorted(:,3).^2) ./ (2.d0.*RSorted(:,1).*RSorted(:,2)) );
  alpha(:,1) = acos( (RSorted(:,2).^2 + RSorted(:,3).^2 - RSorted(:,1).^2) ./ (2.d0.*RSorted(:,2).*RSorted(:,3)) );
  alpha(:,2) = acos( (RSorted(:,1).^2 + RSorted(:,3).^2 - RSorted(:,2).^2) ./ (2.d0.*RSorted(:,1).*RSorted(:,3)) );
  alpha    = alpha .* 180 ./ pi;
%   figure(iFigure)
%   scatter3(RSorted(:,1).*AngToBohr,RSorted(:,2).*AngToBohr,RSorted(:,3).*AngToBohr);
%   iFigure = iFigure + 1;
%   pause
%   figure(iFigure)
%   scatter(RSorted(:,1).*AngToBohr,E_1(:).*KcalMolToEV);
%   iFigure = iFigure + 1;
%   pause
  
  
  ROI = 1.2076;
  AOI = 116.75;
  %FileName1 = strcat('./R',num2str(ROI),'_Alpha',num2str(AOI),'.csv.',num2str(iPES));
  FileName1 = strcat('./Cut_4.csv.',num2str(iPES));
  fileID1   = fopen(FileName1,'w');
  fprintf(fileID1,'R1,R2,R3,RPlot,E\n');
  for i=1:size(RSorted,1)
    j=0;
    if abs(RSorted(i,1) - ROI) < 1.d-3
      if abs(alpha(i,2) - AOI) < 1.d-2
        j=3;
      elseif abs(alpha(i,3) - AOI) < 1.d-2
        j=2;
      end
    elseif abs(RSorted(i,2) - ROI) < 1.d-3
      if abs(alpha(i,1) - AOI) < 1.d-2
        j=3;
      elseif abs(alpha(i,3) - AOI) < 1.d-2
        j=1;
      end
    elseif abs(RSorted(i,3) - ROI) < 1.d-3
      if abs(alpha(i,2) - AOI) < 1.d-2
        j=1;
      elseif abs(alpha(i,1) - AOI) < 1.d-2
        j=2;
      end
    end
    if j~=0
      plot(RSorted(i,j),E_1(i),'ro')
      hold on
      %plot(RSorted(i,j),E_2(i),'ko')
      %hold on
      fprintf(fileID1,'%e,%e,%e,%e,%e\n', RSorted(i,1).*AngToBohr, RSorted(i,2).*AngToBohr, RSorted(i,3).*AngToBohr, RSorted(i,j).*AngToBohr, E_1(i).* KcalMolToEV);
    end
  end
  fclose(fileID1);
  
  RSorted = RSorted .* AngToBohr;
  E_1      = E_1 .* KcalMolToEV;
  E_2      = E_2 .* KcalMolToEV;
  FileName1 = strcat('./RSampled.csv.',num2str(iPES));
  fileID1   = fopen(FileName1,'w');
  %fprintf(fileID1,'R1,R2,R3\n');
  FileName2 = strcat('./ESampled.csv.',num2str(iPES));
  fileID2   = fopen(FileName2,'w');
  %fprintf(fileID2,'E_total,E_CASSCF\n');
  FileName3 = strcat('./RESampled.csv.',num2str(iPES));
  fileID3   = fopen(FileName3,'w');
  fprintf(fileID3,'R1,R2,R3,E\n')
  for i = 1:length(RSorted)
    fprintf(fileID1,'%e,%e,%e\n', RSorted(i,1), RSorted(i,2), RSorted(i,3));
    fprintf(fileID2,'%e,%e\n',    E_1(i),       E_2(i));
    fprintf(fileID3,'%e,%e,%e,%e\n', RSorted(i,1), RSorted(i,2), RSorted(i,3), E_1(i));
  end
  fclose(fileID1);
  fclose(fileID2);
  fclose(fileID3);
  
  clear Col1 Col2 Col3 E_1 E_2 R1 R2 R3 X2 X3 R RSorted alpha
  
end



% for iPES=1:NPESs
% 
%   filename = strcat(WORSPACE_PATH,'/CoarseAIR/run_O3_PlotPES/Test/PlotPES/PES_', num2str(iPES), '/PESFromReadPoints.csv')
%   delimiter = ',';
%   startRow = 2;
%   formatSpec = '%f%f%f%f%[^\n\r]';
%   fileID = fopen(filename,'r');
%   dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%   fclose(fileID);
%   R1 = dataArray{:, 1};
%   R2 = dataArray{:, 2};
%   R3 = dataArray{:, 3};
%   E  = dataArray{:, 4};
%   clearvars filename delimiter startRow formatSpec fileID dataArray ans;
% 
%   filename = strcat(WORSPACE_PATH,'/CoarseAIR/coarseair/dtb/O3/PESs/AbInitioData/E.csv.', num2str(iPES))
%   delimiter = ',';
%   startRow = 2;
%   formatSpec = '%f%f%[^\n\r]';
%   fileID = fopen(filename,'r');
%   dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'HeaderLines' ,startRow-1,  'ReturnOnError', false);
%   fclose(fileID);
%   E_A = dataArray{:, 1};
%   E_B = dataArray{:, 2};
%   clearvars filename delimiter formatSpec fileID dataArray ans;
% 
%   figure(iFigure)
%   scatter(E,E_A)
%   hold on
%   scatter(E,E_B)
%   plot([0,max(E)],[0,max(E)],'-')
%   iFigure = iFigure+1;
%   
%   clear R1 R2 R3 E E_A E_B
%   
% end 