function [Steps]=ComputeMovieSteps(iT, t)    
    
  global NStepsMovie TimeMinMovie TimeMaxMovie

  TimeMovie = logspace(log10(TimeMinMovie(iT)),(log10(TimeMaxMovie(iT))),NStepsMovie);
  iTime = 1;
  for iSteps = 1:NStepsMovie
    while TimeMovie(iSteps) > t(iTime)
      iTime = iTime + 1;
    end
    StepsTemp(iSteps)=iTime;
  end
  
  jSteps   = 2;
  Steps(1) = StepsTemp(1);
  for iSteps = 2:length(StepsTemp)
    if StepsTemp(iSteps-1) ~= StepsTemp(iSteps)
      Steps(jSteps) = StepsTemp(iSteps);
      jSteps = jSteps + 1;
    end
  end

end