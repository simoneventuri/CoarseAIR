close all
clear all
clc

R1Vec  = [[1.2:0.15:1.95],[2.034:0.02:2.114],[2.25:0.10:3.15],3.30,3.60,[4.0:0.75:9.0]];
NR1    = length(R1Vec);

R3Vec  = [[1.2:0.15:2.1],[2.136:0.02:2.216],[2.45:0.10:3.3],3.60,[4.0:0.75:9.0]];
NR3    = length(R3Vec);

Theta2Vec = [10:10:170];
NTheta2   = length(Theta2Vec);

DiatPES  = @(R11,cosTheta2) N2_UMN(R11(1)) + NO_UMN( sqrt( R11(1)^2 +R11(2)^2 - 2.0*R11(1)*R11(2)*cosTheta2 ) ) + NO_UMN(R11(2))
R0  = [2.8,2.8]; 

FileName1 = strcat('./RSampled.csv.1');
fileID1   = fopen(FileName1,'w');
%fprintf(fileID1,'R1,R2,R3\n');
for iTheta2 = 1:NTheta2
  Theta2    = Theta2Vec(iTheta2);
  cosTheta2 = cos(Theta2/180.0*pi);
  for i1 = 1:NR1
    for i3 = 1:NR3
      R1     = R1Vec(i1);
      R3     = R3Vec(i3);
      R2     = sqrt( R1^2 + R3^2 - 2.0*R1*R3*cosTheta2 );
      RR     = [R1,R2,R3];
      if ( ( (R1<(R2+R3)) && (R2<(R1+R3)) && (R3<(R2+R1)) ) )
        fprintf(fileID1,'%e,%e,%e\n', R1, R2, R3 );
      end
    end
  end
end

% for Theta2=[1:1:179]
%   cosTheta2 = cos(Theta2/180.0*pi);
%   R1        = R1Vec(i1);
%   R3        = R3Vec(i3);
%   R2        = sqrt( R1^2 + R3^2 - 2.0*R1*R3*cosTheta2 );
%   %fun       = @(R11,R33)DiatPES(R11,R33,cosTheta2);
%   fun       = @(R11)DiatPES(R11,cosTheta2);
%   RVec      = fminsearch(fun, R0);
%   R1        = RVec(1);
%   R3        = RVec(2);
%   R2        = sqrt( R1^2 + R3^2 - 2.0*R1*R3*cosTheta2 );
%   RVec
%   fprintf(fileID1,'%e,%e,%e\n', R1, R2, R3 );
% end

rlin       = linspace(1.2,7,3000);
[VN2, dV]  = N2_UMN(rlin);
[VN2_, dV] = N2_UMN(R1Vec);
[VNO, dV]  = NO_UMN(rlin);
[VNO_, dV] = NO_UMN(R3Vec);
figure
plot(rlin,VN2)
hold on
plot(R1Vec,VN2_,'o')
plot(rlin,VNO)
plot(R3Vec,VNO_,'o')
