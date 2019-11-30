function [Eint, Rbond, AngMom] = ComputeState(iMol, xx, xxdot)
  
  global Pair_to_Atoms AtomMass
  
  N       = 3;
  NSpace  = 3;
  NEqtVar = 6;
  
  
  for i = 1:2
    x(:,i)    = xx(:,Pair_to_Atoms(iMol,i));
    xdot(:,i) = xxdot(:,Pair_to_Atoms(iMol,i));
  end
  
  xma = AtomMass(Pair_to_Atoms(iMol,1));
  xmb = AtomMass(Pair_to_Atoms(iMol,2));
  
  fa                =   xma / ( xma + xmb );
  fb                =   xmb / ( xma + xmb );
  xcm               =   fa * x(:,1)    + fb * x(:,2);
  xdotcm            =   fa * xdot(:,1) + fb * xdot(:,2);
  itype             =   1;
  
  Eint              =   0.d0;
  Rbond             =   0.d0;
  for i = 1:NSpace
    x(i,:)          =   x(i,:)    - xcm(i);
    xdot(i,:)       =   xdot(i,:) - xdotcm(i);
    Rbond           =   Rbond + ( x(i,1) - x(i,2) )^2;
    Eint            =   Eint  + 5.d-1 * ( xma * xdot(i,1)^2 + xmb * xdot(i,2)^2 );
  end
  Eint;          
  Rbond             =   sqrt(Rbond);
  
  
  AngMomComponents = zeros(3,1);
  for i = 1:2
    iA                  = Pair_to_Atoms(iMol,i); 
    AngMomComponents(1) = AngMomComponents(1) + (x(2,i) * xdot(3,i) - x(3,i) * xdot(2,i)) * AtomMass(iA);
    AngMomComponents(2) = AngMomComponents(2) + (x(3,i) * xdot(1,i) - x(1,i) * xdot(3,i)) * AtomMass(iA);
    AngMomComponents(3) = AngMomComponents(3) + (x(1,i) * xdot(2,i) - x(2,i) * xdot(1,i)) * AtomMass(iA);
  end
  tis           =   sum( AngMomComponents.^2 );
  AngMom        =   sqrt(tis);

end