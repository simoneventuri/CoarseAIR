function [PopVecStoch, PopOverQStoch, PopStoch, NStepsStoch] = ReadPopStoch(iT, QBins)     

    %% (METHOD DEPENDENT)

    global T0_Vec BinnedMolName NBinnedMol NBins DatabasePath OutputPath SystemPath KinMthd iPESStart iPESEnd
    
    for iPES=iPESStart:iPESEnd
    
        for iBinnedMol=1:NBinnedMol

          % Reading Binned Molecules' Bins Populations
          filename = strcat(OutputPath,'/T_',num2str(T0_Vec(iT)),'/output/','/pop_',BinnedMolName(iBinnedMol,:),'.dat.',num2str(iPES))
          delimiter = ' ';
          startRow = 3;
          formatSpec = '%*s%f%*s%*s%[^\n\r]';
          fileID = fopen(filename,'r');
          dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
          fclose(fileID);
          temp=length(dataArray{:, 1});
          PopVecStoch(1:temp,iPES,iBinnedMol) = dataArray{:, 1};
          clearvars delimiter startRow formatSpec fileID dataArray ans;


          iStep = 1;
          iBin  = 1;
          for iVec=1:temp
            if PopVecStoch(iVec,iPES,iBinnedMol) ~= iStep
              PopOverQStoch(iStep,iBin,iPES,iBinnedMol) = PopVecStoch(iVec,iPES,iBinnedMol);
              PopStoch(iStep,iBin,iPES,iBinnedMol)      = PopVecStoch(iVec,iPES,iBinnedMol) .* QBins(iBin,iBinnedMol);
              iBin      = iBin + 1;
            else
              iStep = iStep + 1;
              NBins(iBinnedMol) = iBin-1;
              iBin  = 1;
            end
          end
          NStepsStoch = iStep;

        end
        
    end
    
end