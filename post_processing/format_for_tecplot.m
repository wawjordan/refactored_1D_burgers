function format_for_tecplot(dirname,filename,S)
Nvars = length(S.variables);
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
[mm,nn]=size(DATA);
fprintf(fid,'I=%d\n',mm);
fprintf(fid,zoneStr2);
for j = 1:nn
    for i = 1:mm
        fprintf(fid, [S.dataFmt,'\n'],DATA(i,j));
    end
end
fprintf(fid, '\n');
end
fclose(fid);
end