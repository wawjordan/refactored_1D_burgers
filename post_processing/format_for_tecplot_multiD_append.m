function format_for_tecplot_multiD_append(dirname,filename,S,start)
Nvars = length(S.variables);
dim = S.dim;
if (start)
    headerStr=['TITLE=',S.title,'\n','variables='];
    for i = 1:Nvars-1
        headerStr=[headerStr,'"',char(S.variables{i}),'",'];
    end
    headerStr=[headerStr,'"',char(S.variables{Nvars}),'"\n'];        
end
zoneStr1 = 'ZONE\n';
zoneFmt = ['T = "',S.zoneFmt,'"\n'];
zoneStr2 = ['ZONETYPE=Ordered\n','DT=(',repmat('DOUBLE ',1,Nvars-1),...
    'DOUBLE)\n','DATAPACKING=BLOCK\n'];
name = [dirname,filename];
if (start)
    fid=fopen(name,'wt');
    fprintf(fid,headerStr);
    if (S.aux == true)
        for i = 1:length(S.auxdata)
            fprintf(fid,'DATASETAUXDATA %s = "%s"\n',S.auxdata(i).name,S.auxdata(i).value);
        end
    end
    if (S.cust == true)
        fprintf(fid,'CUSTOMLABELS ');
        for i = 1:length(S.customlabels)-1
            fprintf(fid,'"%s", ',S.customlabels{i});
        end
        fprintf(fid,'"%s"\n',S.customlabels{length(S.customlabels)});
    end
else
    fid=fopen(name,'a');
end

fprintf(fid,zoneStr1);
fprintf(fid,zoneFmt,S.zoneVar);
DATA = S.DATA.dat;
fprintf(fid,zoneStr2);

if length(dim) == 1
    fprintf(fid,'I=%d\n',dim(1));
    fprintf(fid, '\n');
    for VAR = 1:Nvars
        for I = 1:dim(1)
            fprintf(fid, [S.dataFmt,'\n'],DATA(I,VAR));
        end
        fprintf(fid, '\n');
    end
elseif length(dim) == 2
    fprintf(fid,'I=%d\n',dim(1));
    fprintf(fid,'J=%d\n',dim(2));
    for VAR = 1:Nvars
        for J = 1:dim(2)
            for I = 1:dim(1)
                fprintf(fid, [S.dataFmt,'\n'],DATA(I,J,VAR));
            end
        end
        fprintf(fid, '\n');
    end
elseif length(dim) == 3
    fprintf(fid,'I=%d\n',dim(1));
    fprintf(fid,'J=%d\n',dim(2));
    fprintf(fid,'K=%d\n',dim(3));
    for VAR = 1:Nvars
        for K = 1:dim(3)
            for J = 1:dim(2)
                for I = 1:dim(1)
                    fprintf(fid, [S.dataFmt,'\n'],DATA(I,J,K,VAR));
                end
            end
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);
end