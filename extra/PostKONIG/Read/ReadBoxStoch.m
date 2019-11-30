function [t, MolFracs, Ttrans, rho, P, nd, E] = ReadBoxStoch(iT)    

  global T0_Vec NComp OutputPath NPESs
  global PlotMolFracsFlg PlotTemperaturesFlg PlotIntPopFlg PlotKsFlg PlotKsMoleFracsFlg PlotMovieFlg NStepsMovie ...
         TimeMinMovie TimeMaxMovie PlotStepsFlg tInstants PlotStepsSubplotsFlg NSubplots NInstantsPerSubplot iPESStart iPESEnd 


  MolFracs = [];
  Ttrans   = [];
  rho      = [];
  P        = [];
  nd       = [];
  E        = [];
  for iPES=iPESStart:iPESEnd
  
     % Reading Simulation Time
    filename = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/box.dat.', num2str(iPES))
    startRow = 1;
    formatSpec = '%20f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    t = dataArray{:, 1};
    clearvars startRow formatSpec fileID dataArray ans;
    
    % Reading Mole Fractions
    startRow = 1;
    formatSpec = '%*20s';
    for iComp=1:NComp
        formatSpec=strcat(formatSpec,'%20f');
    end
    formatSpec=strcat(formatSpec,'%[^\n\r]');
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    MolFracs(:,1:NComp,iPES) = [dataArray{1:end-1}];

    % Reading Thermodynamic Variables
    startRow = 1;
    formatSpec = '%*20s';
    for iComp=1:NComp
        formatSpec=strcat(formatSpec,'%*20f');
    end
    formatSpec=strcat(formatSpec,'%20f%20f%20f%20f%f%[^\n\r]');
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    Ttrans(:,iPES) = dataArray{:, 1};
    rho(:,iPES)    = dataArray{:, 2};
    P(:,iPES)      = dataArray{:, 3};
    nd(:,iPES)     = dataArray{:, 4};
    E(:,iPES)      = dataArray{:, 5};

  end

end