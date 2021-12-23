clc; clear; close all;
% filename = 'shock_exact_start_alg1_allsave';
filename = 'expand_exact_start';

load(filename);

NX    = OUT.Nx;
N_grids = length(NX);

E_RICH = struct();
for j = 2:N_grids
    Nt = length(OUT.Local_Error_P(j-1).t);
    E_RICH(j).E = cell(Nt,1);
    for i = 1:Nt
        s = 2*(i-1)+1;
         fprintf('%d, %d\n',i,s);
        E_RICH(j).E{i} = -(OUT.Local_Error_P(j).u{s}(1:2:end)-OUT.Local_Error_P(j-1).u{i})/(2^2-1);
    end
end
% for j = 2:N_grids
%     E_RICH(j).E = cell(Nt,1);
%     for i = 1:Nt
%         s1 = 2^(j-2); % 1, 2, 4, ...
%         s2 = 2*s1;    % 2, 4, 8, ...
%         k1 = s1*(i-1)+1; % 1, 2, 3, 4, ...
%         k2 = s2*(i-1)+1; % 1, 3, 5, 7, ...
% %         fprintf('%d, %d, %d %d\n',i,j,k1,k2);
%         E_RICH(j).E{i} = -(OUT.Local_Error_P(j).u{k2}(1:2:end)-OUT.Local_Error_P(j-1).u{k1})/(2^2-1);
%     end
% end
i = 5;
figure(1)
for t = 2:length(E_RICH(i).E)
clf;
plot(OUT.Local_Error_P(i).x,OUT.Local_Error_P(i).E{t},'k')
hold on
plot(OUT.Local_Error_P(i-1).x,E_RICH(i).E{t},'r')
plot(OUT.Local_Error_E(i).x,OUT.Local_Error_E(i).E{t,1},'b')
plot(OUT.Local_Error_E(i).x,OUT.Local_Error_E(i).E{t,2},'g')
legend('Exact','Richardson Extrapolation','ETE','ETE IC')
drawnow;
end
% 
% figure(2)
% plot(OUT.Local_Error_P(i).x,OUT.Local_Error_P(i).Uex{t},'k')
% hold on
% plot(OUT.Local_Error_E(i-1).x,OUT.Local_Error_P(i-1).U{t},'r')
% plot(OUT.Local_Error_E(i).x,OUT.Local_Error_P(i).E{t},'b')
% legend('Exact Solution','Primal (Coarse)','Primal (Fine)')
