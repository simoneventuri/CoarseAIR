clc
iMol = 1;

LevToGroup = zeros(NLevels(iMol),1);

jvToFinaljv = [];
for i = 1:25
  Temp        = i * ones(1,2^min((i-1),3));
  jvToFinaljv = [jvToFinaljv, Temp];
end
% jvToFinaljv = [1:20];
% jvToFinaljv = [jvToFinaljv, 20*ones(1,63)];
jvToFinaljv

for iv = 0:max(Levelvqn(:,iMol))
  for iLevels = 1:NLevels(iMol)
    if Levelvqn(iLevels,iMol) == iv
      jv       = iv;
      ExitFlag = false;
      while LevelEeV(iLevels,iMol) >= vEeVVib(jv+1,iMol) - 1 && ExitFlag == false
        jv = jv + 1;
        if jv >= max(Levelvqn(:,iMol))
          ExitFlag = true;
        end
      end
      LevToGroup(iLevels,1) = iv*(max(jvToFinaljv)) + jvToFinaljv(jv-iv);
      %LevToGroup(iLevels,1) = iv*(max(Levelvqn(:,iMol))+1) + jv;
    end
  end
end

% TempVec    = randperm(max(LevToGroup));
% LevToGroup = TempVec(LevToGroup(:));

FileName = strcat('./TempLevels.csv');
fileID = fopen(FileName,'w');
fprintf(fileID,'id,v,Longitude,Latitude,rIn,EeVVib,EeVRot,Group\n');
for i = 1:NLevels(iMol)
  fprintf(fileID,'%i,%i,%e,%e,%e,%e,%e,%i\n', i, Levelvqn(i), -Leveljqn(i), LevelEeV(i).*10.d0, rIn(i), LevelEeVVib0(i), LevelEeVRot(i), LevToGroup(i));
end
fclose(fileID);

  

A = [LevToGroup, Levelvqn(:,1), Leveljqn(:,1)]