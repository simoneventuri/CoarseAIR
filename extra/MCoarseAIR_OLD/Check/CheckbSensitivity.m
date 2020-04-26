function [iFigure, FinState, Cross, CrossCDF, CrossVar] = CheckbSensitivity(iFigure, iTint, iBins, vqn, jqn)

  global T0_Vec NBins

  filename = strcat('../Test/T_',num2str(T0_Vec(iTint)),'_',num2str(T0_Vec(iTint)),'/Bins_',num2str(iBins),'_0/statistics-bSensitivity.out');
  startRow = 2;
  endRow = 2;
  formatSpec = '%*2s%12s%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  dataArray{1} = strtrim(dataArray{1});
  fclose(fileID);
  TempCell = dataArray{:, 1};
  NRings = str2num(cell2mat(TempCell));
  clearvars startRow endRow formatSpec fileID dataArray ans;
  startRow = 4;
  TempString = '';
  for iRings=1:2*(NRings+1)
    TempString=strcat(TempString,'%18f');
  end
  formatSpec = strcat('%6f%8f%8f%8f%8f%8f',TempString,'%f%[^\n\r]');
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  TempMatrix = [dataArray{1:end-1}];
  clearvars filename startRow formatSpec fileID dataArray ans;

  vqnIni(:) = TempMatrix(:,5);
  jqnIni(:) = TempMatrix(:,4);
  ArrIni(:) = TempMatrix(:,6);
  vqnFin(:) = TempMatrix(:,2);
  jqnFin(:) = TempMatrix(:,1);
  ArrFin(:) = TempMatrix(:,3);
  for iRings=1:NRings
    Cross(:,iRings)    = TempMatrix(:,6+iRings);
    CrossVar(:,iRings) = TempMatrix(:,6+NRings+iRings);
    if iRings > 1
      CrossCDF(:,iRings) = CrossCDF(:,iRings-1) + Cross(:,iRings);
    else
      CrossCDF(:,1) = Cross(:,1) .* 0.d0;
    end
  end
  CrossCDF(:,1:NRings) = CrossCDF(:,1:NRings) ./ CrossCDF(:,end);
  
  filename = '/Users/sventuri/Downloads/elevels_g.dat';
  delimiter = '\t';
  formatSpec = '%*s%*s%*s%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
  fclose(fileID);
  LevelToBin = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;
  LevelToBin = LevelToBin + 1;
  
  for iState = 1:length(vqnFin)
    Found = 0;
    
%     if (ArrFin(iState) == 16) || (ArrFin(iState) == 32) || (ArrFin(iState) == 17) || (ArrFin(iState) == 33)
%       
%       for iLevel = 1:NBins(1)
%         if (vqnFin(iState) == vqn(iLevel,1)) && (jqnFin(iState) == jqn(iLevel,1) )
%           FinState(iState) = iLevel;
%           Found = 1;
%         end
%       end
%       
%     elseif (ArrFin(iState) == 48) || (ArrFin(iState) == 49)
%       
%       for iLevel = 1:NBins(2)
%         if (vqnFin(iState) == vqn(iLevel,2)) && (jqnFin(iState) == jqn(iLevel,2) )
%           FinState(iState) = iLevel + NBins(1);
%           Found = 1;
%         end
%       end
%     
%     end 
     
    if (ArrFin(iState) == 16) || (ArrFin(iState) == 32) || (ArrFin(iState) == 17) || (ArrFin(iState) == 33) || (ArrFin(iState) == 48) || (ArrFin(iState) == 49)
      
      for iLevel = 1:NBins(1)
        if (vqnFin(iState) == vqn(iLevel,1)) && (jqnFin(iState) == jqn(iLevel,1) )
          FinState(iState) = iLevel;
          FinBin(iState)   = LevelToBin(iLevel);
          Found = 1;
        end
      end
    
    end 

    if Found == 0
      FinState(iState) = 0;
      FinBin(iState)   = 0;
    end
    
  end
 
  
  for iRing=1:NRings
    [CrossSort, SortOrder] = sort(Cross(:,iRing),'descend');
    FinBin(SortOrder(:))
    iTemp=1;
    while (iTemp<=length(CrossSort)) && (CrossSort(iTemp) > 5.d-1)
      iTemp=iTemp+1;
    end
    iTemp=iTemp-1;
    figure
    CrossSort(1:iTemp)
    histogram(FinBin(SortOrder(1:min(iTemp,30))))
    LevelToBin(iBins)
    clear CrossSort SortOrder
  end
  
end