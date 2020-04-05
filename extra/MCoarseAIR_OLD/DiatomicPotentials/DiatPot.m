function [Ve, dVe] = DiatPot(r, jqn, iMol)
  
  global AtomMass MoleculesName Pair_to_Atoms
  
 
  if (strcmp(MoleculesName(iMol,:), 'N2_NASA') || strcmp(MoleculesName(iMol,1:2), 'NN_NASA'))
    
    [Vv, dVv] = LeRoy(r);
    
  elseif (strcmp(MoleculesName(iMol,:), 'N2_UMN') || strcmp(MoleculesName(iMol,1:2), 'NN_UMN'))
    
    [Vv, dVv] = N2_UMN(r);
    
  elseif (sum(MoleculesName(iMol,1:2) == 'O2')==2 || sum(MoleculesName(iMol,1:2) == 'OO')==2 )
    
    %[Vv, dVv] = O2_NASA(r);
    [Vv, dVv] = O2_UMN(r);
    
  elseif (sum(MoleculesName(iMol,1:2) == 'CO')==2 || sum(MoleculesName(iMol,1:2) == 'OC')==2 )
    
    [Vv, dVv] = CO_NASA(r);
    
  elseif (sum(MoleculesName(iMol,1:2) == 'CH')==2 || sum(MoleculesName(iMol,1:2) == 'HC')==2 )
    
    [Vv, dVv] = CH_UIUC(r);
    
  elseif (sum(MoleculesName(iMol,1:2) == 'HN')==2 || sum(MoleculesName(iMol,1:2) == 'NH')==2 )
    
    [Vv, dVv] = HN_UIUC(r);
    
  elseif (sum(MoleculesName(iMol,1:2) == 'CN')==2 || sum(MoleculesName(iMol,1:2) == 'NC')==2 )
    
    [Vv, dVv] = CN_UIUC(r);
    
  else
    
    [Vv, dVv] = MorsePot(r);
    
  end
    
  [Vc, dVc]   = CentPot(r, jqn, iMol);
  
  Ve  = (Vv  + Vc * 27.2113839712790);
  %Ve  = (Vv);
  dVe = (dVv + dVc * 27.2113839712790);
  
  % NOTE: REMOVE conversion to eV for Ve;
  % NOTE: DiatPot O2_UMN is in eV;
  
end