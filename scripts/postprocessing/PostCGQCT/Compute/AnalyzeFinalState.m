function [riFin, roFin, jqnFin, vqnFin, iTypeFin] = AnalyzeFinalState(jqnIn, AngMom, Eint, rBond)
   
  global MoleculesName

  Zero  = 0.d0;
  Half  = 5.d-1;
  One   = 1.d0;
  Two   = 2.d0;
  
  rfact         = 0.d0;
  K_c0          = 2.99792458e+08;
  rydberg       = 10973731.568527;
  autime_to_sec = 0.25d0 / ( pi * K_c0 * rydberg );
  TLife         = 100.d0 * 1.d-12 / autime_to_sec;
  epss          = 1.d-4;
  
  jqn0  = AngMom - Half;
  iType = 1;
  ri    = Zero;
  ro    = 100.d0;

  
  [vib, dVv] = LeRoy(rBond);
  Eint       = Eint + vib;
  
  [rMin, VMin, rMax, VMax, VInf] = FindMinMaxDiatPot(jqn0);
  
  if rBond > rMax || Eint < VMin
    
    iType = 3;
    viba  = Half;
    
  else
    
    if Eint < VInf
      iType = 1;      
    else
      iType = 2;
    end
    
    if Eint > VMax 
  
      iType = 3 
      viba  = Half;
      
    else
      
      [rLx, rRx] = TurningPoints( [1.5d0, rMin], [rMin, rMax], Eint, jqn0 );
      plot([rLx, rRx], [Eint, Eint],'-')
      ri = rLx;
      ro = rRx;
      
      [viba] = ActionIntegral( rLx, rRx, Eint, jqn0 );
      
      gam = Zero;
      if Eint > VInf
        
        [rRx, r3] = TurningPoints( [rRx-0.1, rMax], [rMax, 20.d0], Eint, jqn0 );
        %plot([rRx, r3], [Eint, Eint],'-k')
        ri = rRx;
        ro = r3;
        
        Eintp = Eint + epss;
        if ( Eintp > VMax ) 
          Eintp = Eint + Half * ( VMax - Eint );
        end
        dEint        = Eintp - Eint;
        [rLxp, rRxp] = TurningPoints( [1.5d0, rMin], [rMin, rMax], Eintp, jqn0 );
        [Actionp]    = ActionIntegral( rLxp, rRxp, Eintp, jqn0 );
        
        Eintm        = Eint - epss;
        [rLxm, rRxm] = TurningPoints( [1.5d0, rMin], [rMin, rMax], Eintm, jqn0 );
        [Actionm]    = ActionIntegral( rLxm, rRxm, Eintm, jqn0 );
        
        dAdE = Half * ( Actionp - Actionm ) / dEint;
        
        
        [Action] = ActionIntegral( rLx, r3, Eint, jqn0 );
        prob     = exp(-Action);
        gam    = prob / dAdE;
        xlife    = 1 / gam;
        
        if xlife < TLife 
          iType = 4
        end
        
        Action = Action * Half / pi;
        rho2   = gamma(Action) - Action * ( log( abs(Action) ) - One );
        viba   = viba + rho2 * rfact;
      end
      
      viba = viba / (Two * pi);

    end
    
    
  end
  
  iType = iType - 1;
  if iType == 2
    AngMom = Half;
  end
  
  
  jqn = AngMom - Half;
  vqn = viba   - Half;
  jqn = jqn    - Half;
  if mod(floor(jqnIn),2) == 1
    jqn = jqn - One;
    if (round(jqn * Half)  < Zero) && (jqnIn > Zero)
      jqn = jqn + Two;
    end
  end
  jqn = jqn * Half;

  vqn = round(vqn);
  jqn = round(jqn);

  jqn = jqn * Two;
  if mod(floor(jqnIn),2) == 1
    jqn = jqn + One;
  end
  
  
  riFin    = ri;
  roFin    = ro;
  iTypeFin = iType;
  jqnFin   = jqn;
  vqnFin   = vqn;
    
end
