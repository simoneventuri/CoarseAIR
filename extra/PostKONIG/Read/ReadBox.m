function [t, MolFracs, Ttrans, rho, P, nd, E] = ReadBox(iT)    

  global T0_Vec NComp OutputPath
  global PlotMolFracsFlg PlotTemperaturesFlg PlotIntPopFlg PlotKsFlg PlotKsMoleFracsFlg PlotMovieFlg NStepsMovie ...
         TimeMinMovie TimeMaxMovie PlotStepsFlg tInstants PlotStepsSubplotsFlg NSubplots NInstantsPerSubplot

  % Reading Simulation Time
  filename = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/box.dat')
  startRow = 1;
  formatSpec = '%20f%[^\n\r]';
  fileID = fopen(filename,'r');
  dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
  fclose(fileID);
  t = dataArray{:, 1};
  clearvars startRow formatSpec fileID dataArray ans;


  %if PlotMolFracsFlg == 1 || PlotKsMoleFracsFlg == 1 || PlotMovieFlg == 1 || PlotStepsFlg == 1 || PlotStepsSubplotsFlg == 1

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
    MolFracs = [dataArray{1:end-1}];

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
    Ttrans = dataArray{:, 1};
    rho    = dataArray{:, 2};
    P      = dataArray{:, 3};
    nd     = dataArray{:, 4};
    E      = dataArray{:, 5};

  %end

end