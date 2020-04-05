function [Ve] = DiatPotV(r, jqn, iMol)
  
  global AtomMass MoleculesName Pair_to_Atoms
  
  if (strcmp(MoleculesName(iMol,:), 'N2_NASA') || strcmp(MoleculesName(iMol,1:2), 'NN_NASA'))
    
    [Vv, dVv] = LeRoy(r);
    
  elseif (strcmp(MoleculesName(iMol,:), 'N2_UMN') || strcmp(MoleculesName(iMol,1:2), 'NN_UMN'))
    
    [Vv, dVv] = N2_UMN(r);
    
  elseif (strcmp(MoleculesName(iMol,1:2), 'O2') || strcmp(MoleculesName(iMol,1:2), 'OO'))
    
    %[Vv, dVv] = O2_NASA(r);
    [Vv, dVv] = O2_UMN(r);
    
  elseif (strcmp(MoleculesName(iMol,1:2), 'CO') || strcmp(MoleculesName(iMol,1:2), 'OC'))
    
    [Vv, dVv] = CO_NASA(r);
    
  elseif (strcmp(MoleculesName(iMol,1:2), 'CH') || strcmp(MoleculesName(iMol,1:2), 'HC'))
    
    [Vv, dVv] = CH_UIUC(r);
    
  elseif (strcmp(MoleculesName(iMol,1:2), 'HN') || strcmp(MoleculesName(iMol,1:2), 'NH'))
    
    [Vv, dVv] = HN_UIUC(r);
    
  elseif (strcmp(MoleculesName(iMol,1:2), 'CN') || strcmp(MoleculesName(iMol,1:2), 'NC'))
    
    [Vv, dVv] = CN_UIUC(r);
    
  else
    
    [Vv, dVv] = MorsePot(r);
    
  end

  [Vc, dVc]   = CentPot(r, jqn, iMol);
  
  Ve  = (Vv  + Vc * 27.2113839712790);
  %Ve  = (Vv);
  
  % NOTE: REMOVE conversion to eV for Ve;
  % NOTE: DiatPot O2_UMN is in eV;
  
end