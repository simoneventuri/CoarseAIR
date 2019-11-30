function [rLx, rRx, ErrorFlg] = TurningPoints( ExtLx, ExtRx, Eint, jqn, iMol )

  global AtomMass MoleculesName Pair_to_Atoms System
   
  rLx    = fzero(@(x) DiatPotShift(x, jqn, iMol, Eint), ExtLx);
  %VTest1 = DiatPot(rLx, jqn, iMol);
  
  rRx    = fzero(@(x) DiatPotShift(x, jqn, iMol, Eint), ExtRx);
  %VTest2 = DiatPot(rRx, jqn, iMol);
  
  ErrorFlg = false;
    
  tol     = 5.d-13;
  %tol2    = 1.d-15;
%   
%   for iExt = 1:2
% 
%     % BISECTION METHOD
%     if iExt==1
%       x1 = ExtLx(1);
%       x2 = ExtLx(2);
%     else
%       x1 = ExtRx(1);
%       x2 = ExtRx(2);
%     end 
%     [f1, dVe] = DiatPot(x1, jqn, iMol);
%     f1        = Eint - f1;
%     [f2, dVe] = DiatPot(x2, jqn, iMol);
%     f2        = Eint - f2;
%     
%     iter = 0;
%     if f1*f2>0 
%       disp('Bisection Mthd Error!')
%       x1
%       x2
%       f1
%       f2
%       ErrorFlg = 1;
%     else
%       p    = (x1 + x2)/2;
%       [ft, dVe] = DiatPot(p, jqn, iMol);
%       fp = Eint - ft;
%       err = abs(fp);
%       while err > tol
%         iter = iter + 1;
%         if f1*fp < 0  
%           x2 = p;
%         else
%           x1 = p;
%         end
%         p         = (x1 + x2)/2;
%         %pp(iter)  = p;
%         [ft, dVe] = DiatPot(p, jqn, iMol);
%         fp        = Eint - ft;
%         err       = abs(fp);
%       end
%       ErrorFlg = 0;
%     end
%     
%     x  = p;
%     %Rclear pp
%     
%     
% %     % NEWTON METHOD
% %     iTer  = 0;
% %     error = 1.d0;
% %     Vold  = 1.d2;
% %     while abs(error) > tol2
% %       
% %       [Ve, dVe] = DiatPot(x, jqn, iMol);
% %       
% %       V  = Eint - Ve;
% %       
% %       error = V / dVe;
% % 
% %       x = error + x;
% %       
% %       error = V-Vold;
% %       Vold  = V;
% %       iTer = iTer + 1;
% %     end
%  
%     if iExt==1
%       rLx = x;
%     else
%       rRx = x;
%     end
%     
%   end
%   
%   figure
%   rTemp     = linspace(1.5,4,30000);
%   [Ve, dVe] = DiatPot(rTemp', jqn, iMol);
%   plot(rTemp,Ve)
%   hold on
%   plot(rLx, Eint, 'ro')
%   plot(rRx, Eint, 'go')
%   pause
  
end