function [niRatio, ProcessesRatesOverall] = ComputeProcessesRatesOverall(iT, Pop, ProcessesRates, Steps)
    
  global NBins

  iBinnedMol=1;
  size(Pop,1)
  for iStep = 1:size(Pop,1)
    iStep

    SumTemp = sum( Pop(iStep,:,iBinnedMol) );
    
    niRatio(iStep,1:NBins(iBinnedMol),iBinnedMol) = Pop(iStep,1:NBins(iBinnedMol),iBinnedMol) ./ SumTemp;
    
    ProcessesRatesOverall(iStep,1:4,iT) = sum( niRatio(iStep,1:NBins(iBinnedMol),iBinnedMol)' .* max(1.d-100,ProcessesRates(1:NBins(iBinnedMol),1:4,iT)) );

  end

end