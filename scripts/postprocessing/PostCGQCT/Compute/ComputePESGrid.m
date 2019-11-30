function [iPointInPES, PESAngle, RGrid] = ComputePESGrid()

  global AnglesVecGrid RMinGrid RMaxGrid NPointsGrid
  
  RR = linspace(RMinGrid, RMaxGrid, NPointsGrid);
  
  iPoints = 1;
  for Ang=AnglesVecGrid
    jPoints = 1;
    for iR1=1:NPointsGrid
      for iR3=iR1:NPointsGrid
        RGrid(iPoints,1)     = RR(iR1);
        RGrid(iPoints,3)     = RR(iR3);
        RGrid(iPoints,2)     = sqrt( RR(iR1).^2 + RR(iR3).^2 - 2.d0.*RR(iR1).*RR(iR3).*cos(Ang./180.0.*pi) );
        iPointInPES(iPoints) = jPoints;
        PESAngle(iPoints)    = Ang;
        iPoints=iPoints+1;
        jPoints=jPoints+1;
      end
    end
  end
  
end