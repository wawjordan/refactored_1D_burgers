function format_for_tecplot_multiD(dirname,filename,S)
Nvars = length(S.variables);
dim = S.dim;
headerStr=['TITLE=',S.title,'\n','variables='];
for i = 1:Nvars-1
    headerStr=[headerStr,'"',char(S.variables{i}),'",'];
end
headerStr=[headerStr,'"',char(S.variables{Nvars}),'"\n'];

zoneStr1 = 'ZONE\n';
zoneFmt = ['T = "',S.zoneFmt,'"\n'];
zoneStr2 = ['ZONETYPE=Ordered\n','DT=(',repmat('DOUBLE ',1,Nvars-1),...
    'DOUBLE)\n','DATAPACKING=BLOCK\n'];

name = [dirname,filename];
fid=fopen(name,'wt');
fprintf(fid,headerStr);
for k = 1:S.Nzones
    fprintf(fid,zoneStr1);
    fprintf(fid,zoneFmt,S.zoneVar(k));
    DATA = S.DATA(k).dat;
%     [mm,nn]=size(DATA);
    fprintf(fid,zoneStr2);
    if length(dim) == 1
        fprintf(fid,'I=%d\n',dim(1));
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
end
fclose(fid);
end