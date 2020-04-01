function KEq = ComputeThermalRates(NLevels, Levelg, LevelEeV, ProcessesRates)

    global T0_Vec UKb Ue

    iMol = 1;

    ExpVec                       = zeros(size(Levelg));
    ExpVec(1:NLevels(iMol),iMol) = Levelg(1:NLevels(iMol),iMol) .* exp( - LevelEeV(1:NLevels(iMol),iMol) .* Ue ./ (T0_Vec(1) .* UKb) );
    ExpVec(1:NLevels(iMol),iMol) = ExpVec(1:NLevels(iMol),iMol) ./ sum(ExpVec(1:NLevels(iMol),iMol));
%     figure
%     plot(ExpVec(1:NLevels(iMol),iMol))
    
    KEq = zeros(size(ProcessesRates,2),1);
    for iP=1:size(ProcessesRates,2)
        KEq(iP) = sum(ProcessesRates(1:NLevels(iMol),iP) .* ExpVec(1:NLevels(iMol),iMol) );
    end
    
end