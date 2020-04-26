function [RGrid, EGrid, EGridMean, EGridSD] = ComputeStochPES(ReadPESFlg)
  
  global Pair_to_Atoms AtomMass NBins 
  
  global System
  global NHL Network_Folder AbscissaConverter PreLogShift BondOrderFun PIPFun
  global AnglesVecGrid RMinGrid RMaxGrid NPointsGrid
  global iPESStart iPESEnd 
  
  FileNamePES = './StochPES.mat';
  
  NPESTemp = iPESEnd - iPESStart + 1; 
  
  
  if (ReadPESFlg)
    
    load(FileNamePES, 'iPESStart', 'iPESEnd', 'AnglesVecGrid', 'RMinGrid', 'RMaxGrid', 'NPointsGrid', 'iPointInPES', 'PESAngle', 'RGrid', 'EGrid', 'EGridMean', 'EGridSD');
    
  else  

    [iPointInPES, PESAngle, RGrid] = ComputePESGrid();

    PES.Hists.Lambda = [];
    PES.Hists.re     = [];
    PES.Hists.W1     = [];
    PES.Hists.W2     = [];
    PES.Hists.W3     = [];
    PES.Hists.b1     = [];
    PES.Hists.b2     = [];
    PES.Hists.b3     = [];
    PES.Hists.Sigma  = [];
    PES.Hists.Noise  = [];

    [PES.Params.G_MEAN, PES.Params.G_SD] = ReadBNNScales()

    EGrid       = [];
    EGridSum    = [];
    EGridSumSqr = [];
    for iPES=iPESStart:iPESEnd
      iPES
      [PES.Params, PES.Hists] = ReadPESParamsSamples(iPES, PES.Params, PES.Hists);

      EShift          = ComputeBNNPES([100.0, 100.0, 100.0], PreLogShift, PES.Params)
      ETemp(:,1)      = ComputeBNNPES(RGrid, PreLogShift, PES.Params);
      EGrid(:,iPES)   = ETemp - EShift;
      
      
      if iPES==iPESStart
        EGridSum    = EGrid(:,iPES);
        EGridSumSqr = EGrid(:,iPES).^2;
      else
        EGridSum    = EGridSum    + EGrid(:,iPES);
        EGridSumSqr = EGridSumSqr + EGrid(:,iPES).^2;
      end
    end
    EGridMean =      EGridSum    ./ NPESTemp;
    EGridSD   = sqrt(EGridSumSqr ./ NPESTemp - EGridMean.^2);
    
    save(FileNamePES, 'iPESStart', 'iPESEnd', 'AnglesVecGrid', 'RMinGrid', 'RMaxGrid', 'NPointsGrid', 'iPointInPES', 'PESAngle', 'RGrid', 'EGrid', 'EGridMean', 'EGridSD', '-v7.3');
  
  end
  
end