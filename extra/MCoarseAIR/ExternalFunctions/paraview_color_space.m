function [] = paraview_color_space(scheme, FldrName, jsonFlg)

    %colors = eval(scheme);
    %N = length(colors);
    
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
    opts.VariableTypes = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    output = readtable(strcat(FldrName,scheme,'.csv'), opts);
    colors = table2array(output);
    clear opts
    N = size(colors,1);

    if (jsonFlg) 
        
        fid = fopen(strcat(FldrName,scheme,'.json'), 'w');
        fprintf(fid, '[\n');
        fprintf(fid, '    {\n');
        fprintf(fid, '        "IndexedColors" :\n');
        fprintf(fid, '        [\n');
        for i=1:N-1
            for j=1:3
                fprintf(fid, '            %f,\n', colors(i,j));
            end
        end
        fprintf(fid, '            %f,\n', colors(end,1));
        fprintf(fid, '            %f,\n', colors(end,2));
        fprintf(fid, '            %f\n',  colors(end,3));
        fprintf(fid, '        ],\n');
        fprintf(fid, ['        "Name" : "' scheme '",\n']);
        fprintf(fid, '        "NanColor" : \n');
        fprintf(fid, '        [\n');
        fprintf(fid, '            0.80392200000000003,\n');
        fprintf(fid, '            0,\n');
        fprintf(fid, '            0.80392200000000003\n');
        fprintf(fid, '        ]\n');
        fprintf(fid, '    }\n');
        fprintf(fid, ']\n');

        fclose(fid);

        
    else
     
        fid = fopen(strcat(FldrName,scheme,'.xml'), 'w');
        fprintf(fid, '<ColorMaps>\n');
        fprintf(fid, '<ColorMap name="%s" space="RGB">\n', scheme);
        for i=1:N
          x = [(i-1)/(N-1); colors(i,:)'];
          fprintf(fid, '  <Point x="%f" o="1" r="%f" g="%f" b="%f"/>\n', x);
        end
        fprintf(fid, '</ColorMap>\n');
        fprintf(fid, '</ColorMaps>\n');
        fclose(fid);
    
    end

end