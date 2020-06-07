%% O2 Diatomic Potential from Basel (..., 2020)
%
function [V, dV] = O2_Basel(r)

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
    darray1 = [1.0000000000000000        364.17536625912311        5.1196219020132503
               1.1000000000000001       -233.97607003622039        3.5007638586792780
               1.2000000000000000        186.99603321871029        2.3570541559292622
               1.3000000000000000        65.457969582174684        1.5550084704052551
               1.3999999999999999        113.63772881593350       0.99138784336125241
               1.5000000000000000        173.06569791828269       0.59612336885126638
               1.7000000000000000        197.49020582789248       0.13122647858227765
               1.8500000000000001        7.8103646482105704       -4.3510280298733051E-002
               1.9500000000000000       -23.518174159787886      -0.11056706728373911
               2.0499999999999998       -115.86568275189050      -0.15119857062074971
               2.1499999999999999       -142.63049952139795      -0.17320615373574810
               2.2000000000000002       -83.741467105817065      -0.17905724284673852
               2.2500000000000000       -121.88770782990031      -0.18223399504273630
               2.2900000000000000       -98.438846058459532      -0.18316295783472469
               2.3300000000000001       -139.74081802580790      -0.18288057331773189
               2.3799999999999999       -187.80459556193409      -0.18109187058874454
               2.4500000000000002       -393.18208982634650      -0.17642654012274761
               2.5499999999999998       -533.18618662688687      -0.16654651438474843
               2.6499999999999999       -1015.4612425326175      -0.15416007755274563
               2.7999999999999998       -2303.7146152753239      -0.13307582907273741
               3.0000000000000000       -6035.5708423242604      -0.10373041817473450
               3.2500000000000000       -11956.124066045973       -6.9594710854744335E-002
               3.5000000000000000       -22295.493288276270       -4.1696179665734689E-002
               3.7500000000000000        29346.000532199709       -2.2011131192726907E-002
               4.0000000000000000        24832.130491875731       -1.0721461546722821E-002];

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