%% NO Diatomic Potential from Basel (..., 2020)
%
function [V, dV] = NO_Basel(r)

    %%==============================================================================================================
    % 
    % Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
    % 
    % Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
    %
    % Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
    % 
    % This program is free software; you can redistribute it and/or modify it under the terms of the 
    % Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
    % 
    % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
    % without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
    % See the GNU Lesser General Public License for more details. 
    % 
    % You should have received a copy of the GNU Lesser General Public License along with this library; 
    % if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
    % 
    %---------------------------------------------------------------------------------------------------------------
    %%==============================================================================================================
    
    global Param
    
    dk26f1  = 1.0 / 14.0;
    dk26f2  = 1.0 / 18.0;
    
    nda1    = 25;     
    darray1 = [1.0000000000000000        360.15015767133889        4.5159243514525258
               1.1000000000000001       -254.26219186953179        3.0177257098835213
               1.2000000000000000        155.33130631383449        1.9855753035215287
               1.3000000000000000       -8.8268826462189036        1.2755867349725349
               1.3999999999999999        320.24386170362936       0.77788442451051765
               1.5000000000000000        65.770530243624634       0.42308020356753673
               1.7000000000000000        166.66636490753098        9.8374950705135689E-003
               1.8500000000000001       -13.207349385241425      -0.13818715217647082
               1.9500000000000000       -52.709058443412957      -0.19078931138446364
               2.0499999999999998       -130.48313321941350      -0.21888487313847804
               2.1000000000000001       -79.329839945469573      -0.22613555041948530
               2.1499999999999999       -131.71485979566003      -0.22985823254148841
               2.1829999999999998       -16.735953013523069      -0.23068535972748805
               2.2000000000000002       -164.56214560535074      -0.23066827985846317
               2.2500000000000000       -134.38629628245667      -0.22909021963948817
               2.2999999999999998       -412.13044407117184      -0.22557064247448011
               2.3999999999999999       -563.18682481561325      -0.21416819001248655
               2.5000000000000000       -1044.8415846964670      -0.19885616194449085
               2.6499999999999999       -1760.0656149280139      -0.17215852991546399
               2.7999999999999998       -3059.4999610880318      -0.14421321594448955
               3.0000000000000000       -7576.5110352881920      -0.10880003059747878
               3.2500000000000000       -8012.9364314745035       -7.1497695508469405E-002
               3.5000000000000000        2579.2270840447754       -4.4268030907488765E-002
               3.7500000000000000        8889.7347039619199       -2.6635780788467400E-002
               4.0000000000000000        10454.776546313966       -1.6095593593490776E-002];

    ener = r .* 0.0; 
    der  = r .* 0.0;
    for kk = 1:nda1
        xi = darray1(kk,1); 
        
        for iR = 1:length(r)
            xl       = 0.0;
            xs       = 0.0;
            ddrker26 = 0.0;
            drker26  = 0.0;
            
            x = r(iR);
            if (x < xi)
                xl = xi;
                xs = x;
                ddrker26 = -dk26f2 ./ xi.^8;
            else
                xl = x;
                xs = xi;
                ddrker26 = (-7.0 .* dk26f1 ./ x.^8) + (8.0 .* dk26f2 .* xi  ./ x.^9);
            end
            drker26      = (       dk26f1 ./ xl.^7) - (       dk26f2 .* xs ./ xl.^8);

            ener(iR) = ener(iR) +  drker26 * darray1(kk,2);
            der(iR)  = der(iR)  + ddrker26 * darray1(kk,2);
            
        end
    end
    
    V  = ener .* Param.EhToeV;
    dV = der  .* Param.EhToeV;
    
end