function [G_MEAN, G_SD] = ReadBNNScales()

  global Network_Folder
  
  filename = strcat(Network_Folder,'/ScalingValues.csv')
  delimiter = ',';
  startRow = 2;
  endRow = 3;
  formatSpec = '%f%f%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  ScalingValues = [dataArray{1:end-1}];
  clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;
  
  G_MEAN = ScalingValues(1,:);
  G_SD   = ScalingValues(2,:);

end