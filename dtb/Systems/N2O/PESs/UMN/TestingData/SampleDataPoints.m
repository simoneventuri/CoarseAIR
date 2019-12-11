close all
clear all
clc

R1Vec  = [[1.2:0.15:2.05],2.074,[2.1:0.15:3],[3.5:0.5:6.0]];
NR1    = length(R1Vec);

R3Vec  = [[1.2:0.15:2.05],2.176,[2.1:0.15:3],[3.5:1.2:9.5]];
NR3    = length(R3Vec);

Theta2Vec = [10:10:170];
NTheta2   = length(Theta2Vec);


FileName1 = strcat('./RSampled.csv.1');
fileID1   = fopen(FileName1,'w');
%fprintf(fileID1,'R1,R2,R3\n');
for i1 = 1:NR1
  for i3 = 1:NR3
    for iTheta2 = 1:NTheta2
      R1     = R1Vec(i1);
      R3     = R3Vec(i3);
      Theta2 = Theta2Vec(iTheta2);
      R2     = sqrt( R1^2 + R3^2 - 2.0*R1*R3*cos(Theta2/180.0*pi) );
      RR     = [R1,R2,R3];
      if ( ( (R1<(R2+R3)) && (R2<(R1+R3)) && (R3<(R2+R1)) ) )
        fprintf(fileID1,'%e,%e,%e\n', R1, R2, R3 );
      end
    end
  end
end

rlin      = linspace(1.2,7,3000);
[VN2, dV] = N2_UMN(rlin);
[VNO, dV] = NO_UMN(rlin);
figure
plot(rlin,VN2)
hold on
plot(rlin,VNO)
