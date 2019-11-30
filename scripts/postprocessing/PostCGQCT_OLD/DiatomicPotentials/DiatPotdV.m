function [dVe] = DiatdVPot(r, jqn, iMol)
  
  global AtomMass MoleculesName Pair_to_Atoms
  
  if (MoleculesName(iMol,1:2) == 'N2' || MoleculesName(iMol,1:2) == 'NN')
    
    [Vv, dVv] = LeRoy(r);
    
  elseif (MoleculesName(iMol,1:2) == 'O2' || MoleculesName(iMol,1:2) == 'OO')
    
    %[Vv, dVv] = O2_NASA(r);
    [Vv, dVv] = O2_UMN(r);
    
  elseif (MoleculesName(iMol,1:2) == 'CO' || MoleculesName(iMol,1:2) == 'OC')
    
    [Vv, dVv] = CO_NASA(r);
    
  elseif (MoleculesName(iMol,1:2) == 'CH' || MoleculesName(iMol,1:2) == 'HC')
    
    [Vv, dVv] = CH_UIUC(r);
    
  elseif (MoleculesName(iMol,1:2) == 'HN' || MoleculesName(iMol,1:2) == 'NH')
    
    [Vv, dVv] = HN_UIUC(r);
    
  elseif (MoleculesName(iMol,1:2) == 'CN' || MoleculesName(iMol,1:2) == 'NC')
    
    [Vv, dVv] = CN_UIUC(r);
    
  end

  [Vc, dVc]   = CentPot(r, jqn, iMol);

  dVe = (dVv + dVc * 27.2113839712790);
  
  % NOTE: REMOVE conversion to eV for Ve;
  % NOTE: DiatPot O2_UMN is in eV;
  
end