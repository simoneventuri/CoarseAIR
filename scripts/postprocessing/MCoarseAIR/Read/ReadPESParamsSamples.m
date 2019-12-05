function [Params, Hist] = ReadPESParamsSamples(iSample, Params, Hist)
  
  global NHL Network_Folder
  
  FolderName  = 'CalibratedParams';
  FolderName2 = strcat(Network_Folder,'/',FolderName);

  
  filename = strcat(FolderName2,'/Lambda.csv.',num2str(iSample))
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  LambdaTemp = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;
  Params.Lambda = ones(3,1).*LambdaTemp;
  
  filename = strcat(FolderName2,'/re.csv.',num2str(iSample));
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  reTemp = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;
  Params.re = ones(3,1).*reTemp;
  
  
  filename = strcat(FolderName2,'/b1.csv.',num2str(iSample));
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.b1 = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;

  filename = strcat(FolderName2,'/W1.csv.',num2str(iSample));
  delimiter = ',';
  formatSpec = '';
  for i = 1:NHL(2)
    formatSpec = strcat(formatSpec,'%f');
  end
  formatSpec = strcat(formatSpec,'%[^\n\r]');
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.W1 = [dataArray{1:end-1}];
  clearvars filename delimiter formatSpec fileID dataArray ans;


  filename = strcat(FolderName2,'/b2.csv.',num2str(iSample));
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.b2 = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;

  filename = strcat(FolderName2,'/W2.csv.',num2str(iSample));
  delimiter = ',';
  formatSpec = '';
  for i = 1:NHL(3)
    formatSpec = strcat(formatSpec,'%f');
  end
  formatSpec = strcat(formatSpec,'%[^\n\r]');
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.W2 = [dataArray{1:end-1}];
  clearvars filename delimiter formatSpec fileID dataArray ans;


  filename = strcat(FolderName2,'/b3.csv.',num2str(iSample));
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.b3 = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;

  filename = strcat(FolderName2,'/W3.csv.',num2str(iSample));
  delimiter = ',';
  formatSpec = '';
  for i = 1:NHL(4)
    formatSpec = strcat(formatSpec,'%f');
  end
  formatSpec = strcat(formatSpec,'%[^\n\r]');
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.W3 = [dataArray{1:end-1}];
  clearvars filename delimiter formatSpec fileID dataArray ans;


  filename = strcat(FolderName2,'/Sigma.csv.',num2str(iSample));
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.Sigma = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;
  
  
  filename = strcat(FolderName2,'/Noise.csv.',num2str(iSample));
  delimiter = '';
  formatSpec = '%f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
  fclose(fileID);
  Params.Noise = dataArray{:, 1};
  clearvars filename delimiter formatSpec fileID dataArray ans;
  
  
  Hist.Lambda(iSample,:) = Params.Lambda;
  Hist.re(iSample,:)     = Params.re;
  
  Hist.W1(iSample,:,:) = Params.W1;
  Hist.W2(iSample,:,:) = Params.W2;
  Hist.W3(iSample,:,:) = Params.W3;
  
  Hist.b1(iSample,:)   = Params.b1;
  Hist.b2(iSample,:)   = Params.b2;
  Hist.b3(iSample,:)   = Params.b3;
  
  Hist.Sigma(iSample)  = Params.Sigma;
  Hist.Noise(iSample)  = Params.Noise;
  
end