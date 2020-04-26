function [StpInstants]=ComputeStepInstants(t)    

    global tInstants

    StpInstants = ones(1,length(tInstants));
    StpInstants
    for iInstants = 1:length(tInstants)
        while t(StpInstants(iInstants)) < tInstants(iInstants)
            StpInstants(iInstants) = StpInstants(iInstants) + 1;
        end
        StpInstants(iInstants+1) = StpInstants(iInstants);
    end

end