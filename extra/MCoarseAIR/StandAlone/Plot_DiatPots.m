close all
clear all
clc

global Input Syst Param


Input.WORKSPACE_PATH = '/home/venturi/WORKSPACE/'
Input.DATA_PATH      = '/Spebus/PESConstruction/AbInitio_Data/CNH/';%'/Spebus/PESConstruction/AbInitio_Data/O3/UMN_AbInitio/';

Input.SystNameLong = 'O3_UMN';
%Input.SystNameLong = 'N4_NASA';
%Input.SystNameLong   = 'CNH_UIUC';
Input.FigureFormat   = 'PrePrint';
Input.iFig           = 101;

Input.RConv          = 1.0;
Input.EConv          = 1.0;
Input.EShift         = 0.0;%120.311d0



Syst.NameLong = Input.SystNameLong;
Syst          = Initialize_ChemicalSyst(Syst)
Initialize_Parameters()


figure(Input.iFig)
RVec = linspace(1e-10, 8.0, 3000);
for iP=1:3
    [Ve, dVe] = DiatPot(RVec, 0, 1);
    plot(RVec, Ve)
    hold on
    clear Ve dVe
    ylim([-10,10])
end


aa = [2.000000,3.477032,1.477032];

clear( Syst.Molecule(iP).DiatPot )   
[V1, dVv] = feval(Syst.Molecule(1).DiatPot, aa(1));
clear( Syst.Molecule(iP).DiatPot )   
[V2, dVv] = feval(Syst.Molecule(2).DiatPot, aa(2));
clear( Syst.Molecule(iP).DiatPot )   
[V3, dVv] = feval(Syst.Molecule(3).DiatPot, aa(3));
VV = V1+V2+V3

RR = linspace(1.2,10.0,150);
for iAng = [0:5:180]
    
    File   = strcat(Input.WORKSPACE_PATH, Input.DATA_PATH, '/PES_1/VDiat.csv.', num2str(iAng));
    fileID = fopen(File,'w');
    fprintf(fileID,'R1,R2,R3,E\n');
    
    for R1 = RR
        for R3 = RR
            
            R2    = sqrt(R1.^2 + R3.^2 - 2.d0.*R1.*R3.*cos(iAng ./ 180.0 .* pi));
            R     = [R1,R2,R3];
            VDiat = 0.0;
            for iP=1:3
                clear( Syst.Molecule(iP).DiatPot )   
                [Ve, dVv] = feval(Syst.Molecule(iP).DiatPot, R(iP));
                VDiat     = VDiat + Ve;
            end

            fprintf(fileID,'%f,%f,%f,%f\n', R1, R2, R3, VDiat);
            
        end
    end
    
    fclose(fileID);
    
end