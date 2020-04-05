function Write3DLinesToVTK(dir, file, X, Y, Z, ColorVec, cmap)

  NLines  = size(X,1);
  NPoints = NLines .* 2;
  NColors = size(cmap,1);
   
  fid=fopen([dir 'UniformBinning' strcat('.vtk')],'wt');
  fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n');

  fprintf(fid,'POINTS %d float\n',NPoints);
  for i = 1:NLines
      fprintf(fid,'%3.7f %3.7f %3.7f\n',[ 0.d0 0.d0 Z(i,1)]);
      fprintf(fid,'%3.7f %3.7f %3.7f\n',[10.d0 0.d0 Z(i,2)]);
  end

  fprintf(fid,'\nLINES %d %d\n',NPoints,3*NPoints);
  for i = 1:NLines
    fprintf(fid,'2 %d %d\n',(i-1)*2,(i-1)*2+1);
    fprintf(fid,'2 %d %d\n',(i-1)*2+1,(i-1)*2);
  end
   
  fprintf(fid,'\nPOINT_DATA %d \n',NPoints);
  fprintf(fid,'SCALARS element float 2\n');
  fprintf(fid,'LOOKUP_TABLE mytable \n');
  for i = 1:NPoints
    iLine = floor((i-1)/2);
    fprintf(fid,'%3.1f\n',ColorVec(iLine+1));%./max(ColorVec));
    fprintf(fid,'%3.1f\n',ColorVec(iLine+1));%./max(ColorVec));
  end
  
  fprintf(fid,'LOOKUP_TABLE my_table %d \n',NColors);
  for i = 1:NColors
    fprintf(fid,'%3.3f %3.3f %3.3f\n',cmap(i,:));
  end
  
  
  fid=fopen([dir file strcat('.vtk')],'wt');
  fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n');

  fprintf(fid,'POINTS %d float\n',NPoints);
  for i = 1:NLines
      fprintf(fid,'%3.7f %3.7f %3.7f\n',[X(i,1) Y(i,1) Z(i,1)]);
      fprintf(fid,'%3.7f %3.7f %3.7f\n',[X(i,2) Y(i,2) Z(i,2)]);
      
      Mat((i-1)*2+1,1) = X(i,1);
      Mat((i-1)*2+2,1) = X(i,2);
      Mat((i-1)*2+1,2) = Y(i,1);
      Mat((i-1)*2+2,2) = Y(i,2);
      Mat((i-1)*2+1,3) = Z(i,1);
      Mat((i-1)*2+2,3) = Z(i,2);
      
  end

  fprintf(fid,'\nLINES %d %d\n',NPoints,3*NPoints);
  for i = 1:NLines
    fprintf(fid,'2 %d %d\n',(i-1)*2,(i-1)*2+1);
    fprintf(fid,'2 %d %d\n',(i-1)*2+1,(i-1)*2);
  end
   
  fprintf(fid,'\nPOINT_DATA %d \n',NPoints);
  fprintf(fid,'SCALARS element float 2\n');
  fprintf(fid,'LOOKUP_TABLE mytable \n');
  for i = 1:NPoints
    iLine = floor((i-1)/2);
    fprintf(fid,'%3.1f\n',ColorVec(iLine+1));%./max(ColorVec));
    fprintf(fid,'%3.1f\n',ColorVec(iLine+1));%./max(ColorVec));
  end
  
  fprintf(fid,'LOOKUP_TABLE my_table %d \n',NColors);
  for i = 1:NColors
    fprintf(fid,'%3.3f %3.3f %3.3f\n',cmap(i,:));
  end
  
  
  FileName = strcat(dir,'DiatPot.csv');
  csvwrite(FileName,Mat);
  
end