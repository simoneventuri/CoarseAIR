function [niRatioStoch, ProcessesRatesOverallStoch] = ComputeProcessesRatesOverallStoch(iT, PopStoch, ProcessesRatesStoch)
    
  global NBins iPESStart iPESEnd

  iBinnedMol=1;
  size(PopStoch,1)
  for iStep = 1:size(PopStoch,1)
    
    for iPES=iPESStart:iPESEnd

        SumTemp = sum( PopStoch(iStep,:,iPES,iBinnedMol) );

        niRatioStoch(iStep,1:NBins(iBinnedMol),iPES,iBinnedMol) = PopStoch(iStep,1:NBins(iBinnedMol),iPES,iBinnedMol) ./ SumTemp;

        ProcessesRatesOverallStoch(iStep,1:4,iPES) = sum( niRatioStoch(iStep,1:NBins(iBinnedMol),iPES,iBinnedMol)' .* max(1.d-100,ProcessesRatesStoch(1:NBins(iBinnedMol),1:4,iPES)) );
        
    end

  end

end