close all
clear all
clc

R1Vec  = [[1.2:0.15:2.05],2.074,[2.1:0.15:3],[3.5:0.5:6.0]];
NR1    = length(R1Vec);

R3Vec  = linspace(1.5,12,20);
NR3    = length(R3Vec);

Theta2Vec = [10:15:170];
NTheta2   = length(Theta2Vec);


FileName1 = strcat('./N3_Grid.csv');
fileID1   = fopen(FileName1,'w');
fprintf(fileID1,'R1,R2,R3\n');
for i1 = 1:NR1
  for i3 = 1:NR3
    for Theta2 = 1:NTheta2
      R1 = R1Vec(i1);
      R3 = R3Vec(i3);
      R2 = sqrt( R1^2 + R3^2 + 2.0*R1*R3*cos(Theta2/180.0*pi) );
      if ( (R1<(R2+R3)) && (R2<(R1+R3)) && (R3<(R2+R1)) )
        fprintf(fileID1,'%e,%e,%e\n', R1, R2, R3 );
      end
    end
  end
end

rlin    = linspace(1.2,7,3000);
[V, dV] = LeRoy(rlin);
figure
plot(rlin,V)