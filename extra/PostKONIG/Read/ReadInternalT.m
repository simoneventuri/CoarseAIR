function [Tint] = ReadInternalT(iT)    
    
    global T0_Vec NBinnedMol OutputPath
    
    
    filename = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/Tint.dat');        startRow = 3;
    formatSpec = '%*20s';
    for iBinnedMol=1:NBinnedMol
        formatSpec=strcat(formatSpec,'%20f');
    end
    formatSpec=strcat(formatSpec,'%[^\n\r]');
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    Tint = [dataArray{1:end-1}];
    clearvars filename startRow formatSpec fileID dataArray ans;

end