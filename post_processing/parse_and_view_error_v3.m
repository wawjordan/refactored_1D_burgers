%% Parse & View Error V2 (02/11/2022)
clc; clear; close all;

norm = 1;

src_dirname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\new\results',...
    '\02_11\'];
src_fnames = dir([src_dirname,'unsteady_shock_GS*.mat']);
N = length(src_fnames);

Gendata = struct();

Exdata = struct();
Pxdata = struct();

Egdata = struct();
Pgdata = struct();

for i = 1:N
    [~,~,S] = parse_IC_OOA_new([src_dirname,src_fnames(i).name]);
    Gendata(i).N = S.nNodes;
    Gendata(i).t = S.t;
    
    Exdata(i).PRI = S.PRI_Ex(:,:,norm);
    Exdata(i).ETE = S.ETE_Ex(:,:,norm);
    Exdata(i).IC = S.IC_Ex(:,:,norm);

    Pxdata(i).PRI = S.PRI_px(:,:,norm);
    Pxdata(i).ETE = S.ETE_px(:,:,norm);
    Pxdata(i).IC = S.IC_px(:,:,norm);
    
    Egdata(i).PRI = S.PRI_Eg(:,norm);
    Egdata(i).ETE = S.ETE_Eg(:,norm);
    Egdata(i).IC = S.IC_Eg(:,norm);
    
    Pgdata(i).PRI = S.PRI_pg(:,norm);
    Pgdata(i).ETE = S.ETE_pg(:,norm);
    Pgdata(i).IC = S.IC_pg(:,norm);
end

% Plotting

% figure(1);
% subplot(1,2,1)
% hold on
% for i = 1:M2
%     for k = 1:O2
%         plot(dx,PRI_Ex(:,times(i),norms(k)),'k');
%         plot(dx,ETE_Ex(:,times(i),norms(k)),'b');
%         plot(dx, IC_Ex(:,times(i),norms(k)),'r');
%     end
% end
% for k = 1:O2
%     plot(dx,PRI_Eg(:,norms(k)),'m');
%     plot(dx,ETE_Eg(:,norms(k)),'g');
%     plot(dx, IC_Eg(:,norms(k)),'c');
% end
% hold off
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% 
% subplot(1,2,2)
% hold on
% for i = 1:M2
%     for k = 1:O2
%         plot(dx(2:end),PRI_px(2:end,times(i),norms(k)),'k');
%         plot(dx(2:end),ETE_px(2:end,times(i),norms(k)),'b');
%         plot(dx(2:end), IC_px(2:end,times(i),norms(k)),'r');
%     end
% end
% for k = 1:O2
%     plot(dx(2:end),PRI_pg(2:end,norms(k)),'m');
%     plot(dx(2:end),ETE_pg(2:end,norms(k)),'g');
%     plot(dx(2:end), IC_pg(2:end,norms(k)),'c');
% end
% hold off
% set(gca,'xscale','log')