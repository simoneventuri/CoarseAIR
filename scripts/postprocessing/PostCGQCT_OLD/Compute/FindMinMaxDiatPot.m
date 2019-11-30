function [rMin, VMin, rMax, VMax, VInf] = FindMinMaxDiatPot(jqn0, iMol, rMinOld, rMaxOld)

  global MoleculesName AtomMass Pair_to_Atoms
  
  rPlotMin = 1.5d0;
  rPlotMax = 50.d0;
  [VInf, dVTemp] = DiatPot(rPlotMax, jqn0, iMol); 
  rMin = fzero(@(x) DiatPotdV(x, jqn0, iMol), rMinOld);
  VMin = DiatPot(rMin, jqn0, iMol);
  
  rMax = fzero(@(x) DiatPotdV(x, jqn0, iMol), rMaxOld);
  VMax = DiatPot(rMax, jqn0, iMol);
  
  
% %   
% %   rNPoints = 1000;
% %   
% %   
% %   rTemp  = rPlotMin; 
% %   dVTemp = -1.d0;
% %   while dVTemp < 0
% %     rTempOld = rTemp;
% %     rTemp    = rTemp + 0.01d0;
% %     [VTemp, dVTemp] = DiatPot(rTemp, jqn0, iMol); 
% %   end
% %   
% %   rvec = linspace(rTempOld, rTemp, 3000);
% %   for ir = 1:3000
% %     [Ve(ir), dVe(ir)] = DiatPot(rvec(ir), jqn0, iMol);
% %   end
% %   
% %   [VMin, MinIndx] = min(Ve);
% %   rMin = rvec(MinIndx);  
% %  
% %   
% %   rTemp  = rMin; 
% %   dVTemp = 1.d0;
% %   while dVTemp > 0 && rTemp < rPlotMax
% %     rTempOld = rTemp;
% %     rTemp    = rTemp + 0.0001d0;
% %     [VTemp, dVTemp] = x; 
% %   end
% %   
% %   rvec = linspace(rTempOld, rTemp, 10000);
% %   for ir = 1:10000
% %     [Ve(ir), dVe(ir)] = DiatPot(rvec(ir), jqn0, iMol);
% %   end
% %   %figure
% %   %plot(rvec,Ve)
% %   %pause
% %   
% %   [VMax, MaxIndx] = max(Ve);
% %   rMax = rvec(MaxIndx);  
% %  
  
%   clear rvec Ve
%   rvec = linspace(rPlotMin, rPlotMax, rNPoints);
%   for ir = 1:rNPoints
%     [Ve(ir), dVe(ir)] = DiatPot(rvec(ir), jqn0, iMol);
%   end
%   plot(rvec,Ve)
%   hold on
%   plot(rMin,VMin,'ro')
%   plot(rMax,VMax,'ro')
   
end