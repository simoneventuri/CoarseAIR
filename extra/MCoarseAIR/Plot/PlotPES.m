function [iFigure] = PlotPES(iFigure)
  
  global Angles MinPES MaxPES Pair_To_BinnedMolecule MoleculesName

  for iA = 1:length(Angles)
    
    filename = strcat('../Test/PlotPES/Info.dat');

    startRow = 2;
    endRow = 2;
    formatSpec = '%20f%20f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    NPoints(1) = dataArray{:, 1};
    NPoints(2) = dataArray{:, 2};
    clearvars startRow endRow formatSpec dataArray ans;

    startRow = 4;
    endRow = 4;
    fileID = fopen(filename,'r');
    formatSpec = '%20f%20f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    RMin(1) = dataArray{:, 1};
    RMin(2) = dataArray{:, 2};
    clearvars startRow endRow formatSpec dataArray ans;

    startRow = 6;
    endRow = 6;
    fileID = fopen(filename,'r');
    formatSpec = '%20f%20f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    RMax(1) = dataArray{:, 1};
    RMax(2) = dataArray{:, 2};
    clearvars startRow endRow formatSpec fileID dataArray ans;
    clearvars filename startRow formatSpec fileID dataArray ans;

    filename = strcat('../Test/PlotPES/',num2str(Angles(iA)),'Degree_3.dat');
    startRow = 2;
    formatSpec = '%20f%20f%20f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    R1 = dataArray{:, 1};
    R2 = dataArray{:, 2};
    R3 = dataArray{:, 3};
    V = dataArray{:, 4};
    clearvars filename startRow formatSpec fileID dataArray ans;

    x = linspace(RMin(1),RMax(2),NPoints(1));
    y = linspace(RMin(2),RMax(2),NPoints(2));
    [X,Y] = meshgrid(x,y);
    iTot=1;
    for i=1:NPoints(1)
      for j=1:NPoints(2)
        Z(i,j) = V(iTot);
        iTot=iTot+1;
      end
    end
    
    
    figure(iFigure) 
    surf(X,Y,Z)
    xlim([MinPES(1), MaxPES(1)]);
    ylim([MinPES(2), MaxPES(2)]);
    xstring = strcat('$R_{',MoleculesName(Pair_To_BinnedMolecule(1)),'}',{' '},'[a_{0}]$');
    ystring = strcat('$R_{',MoleculesName(Pair_To_BinnedMolecule(3)),'}',{' '},'[a_{0}]$');
    xlabel(xstring);
    ylabel(ystring);
    zlabel(['$Energy [Eh]$']);
    str_title = strcat('Angle = ',{' '},num2str(Angles(iA)));
    title(str_title);
    set(gca,'FontSize',20, 'FontName','Palatino','TickDir','out','TickLabelInterpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');
    az = 130;
    el = 15;
    view(az, el);
    iFigure = iFigure + 1;
    
  end

end