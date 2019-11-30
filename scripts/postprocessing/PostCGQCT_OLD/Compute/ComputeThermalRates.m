function KEq = ComputeThermalRates(Levelg, LevelEeV, Ue, ProcessesRates)

    global NLevels T0_Vec UKb Ue

    iMol = 1;

    ExpVec(1:NLevels(iMol),iMol)  = Levelg(1:NLevels(iMol),iMol) .* exp( - LevelEeV(1:NLevels(iMol),iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
    ExpVec(:,iMol) = ExpVec(:,1) ./ sum(ExpVec(:,iMol));

    for iP=1:size(ProcessesRates,2)
        KEq(iP) = sum(ProcessesRates(1:NLevels(iMol),iP) .* ExpVec(1:NLevels(iMol),iMol) );
    end
    
end