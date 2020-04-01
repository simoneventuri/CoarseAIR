function Write3DLine2VTK(dir, file, iFile, path3D, iColor)
%
% This function takes as input a 3D line and export
% it to an ASCII VTK file which can be oppened with the viewer Paraview.
%
% Input :
%           "file" is the name without extension of the file (string).
%           "path3D" is a list of 3D points (nx3 matrix)
%           "dir" is the path of the directory where the file is saved (string). (Optional)
%
% Simple example :
%
%   Xi=linspace(0,10,200)';
%   Yi=[zeros(100,1);linspace(0,20,100)'];
%   Zi=5*sin(5*Xi);
%   export3Dline2VTK('ExampleLine',[Xi Yi Zi])
%
% David Gingras, January 2009

if nargin==2
    dir=cd;  
end

if ~strcmpi(dir(end),'/'); dir(end+1)='/'; end

X=path3D(:,1);
Y=path3D(:,2);
Z=path3D(:,3);
nbWaypoint=length(X);
if mod(nbWaypoint,3)==1
    X(nbWaypoint+1:nbWaypoint+2,1)=[0;0];
    Y(nbWaypoint+1:nbWaypoint+2,1)=[0;0];
    Z(nbWaypoint+1:nbWaypoint+2,1)=[0;0];
elseif mod(nbWaypoint,3)==2
    X(nbWaypoint+1,1)=0;
    Y(nbWaypoint+1,1)=0;
    Z(nbWaypoint+1,1)=0;
end
nbpoint=length(X);

fid=fopen([dir file strcat('.vtk.',num2str(iFile))],'wt');
fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n');
fprintf(fid,'POINTS %d float\n',nbpoint);

for i=length(X)
  if mod(jqnVec(i),2)==0  
    fprintf(fid,'%3.7f %3.7f %3.7f\n',[X(i,1) Y(i,1) Z(i,1)]');
  end
end

if mod(nbWaypoint,2)==0
    nbLine=2*nbWaypoint-2;
else
    nbLine=2*(nbWaypoint-1);
end

ass1=zeros(nbLine,1);
ass2=zeros(nbLine,1);

ass1(1)=0;
ass2(1)=1;
ass1(end)=nbLine/2;
ass2(end)=nbLine/2-1;

ass1(2:2:nbLine-1)=1:nbLine/2-1;
ass1(3:2:nbLine)=1:nbLine/2-1;
ass2(2:2:nbLine-1)=ass1(2:2:nbLine-1)-1;
ass2(3:2:nbLine)=ass1(3:2:nbLine)+1;

fprintf(fid,'\nLINES %d %d\n',nbLine,3*nbLine);
fprintf(fid,'2 %d %d\n',ass1',ass2');

fclose(fid);