%% Proyecto 6 - Derivadas topológicas
clear all; clc; close all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
clear red blue green purple yellow orange
% ------------------------------------
caso='3';
caso2='1';
fpath = strcat('D:\Google Drive\M2i\Problemas inversos\Proyecto_6\Latex\Figures');
filename = 'Figure';
stringreceptores = ' 1:20';
N=100;
DT_prom = zeros(N,N);
for t=1:11
    load(strcat('DT',num2str(t),'.mat'));
    DT_prom = DT_prom + MatDT/abs((min(min(MatDT))));
end
DT_prom = DT_prom./11;
%%
[X,Y]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));

fig = figure;
set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
box on;hold on; axis tight;
% set(gca,'xTick',[ ],'yTick',[ ])
xlabel('$x$ [-]','FontSize',20,'interpreter','latex');
ylabel('$y$ [-]','FontSize',20,'interpreter','latex');
tit = strcat('Receptores:',' ',stringreceptores);
title(tit,'interpreter','latex','FontSize',20);
contourf(X,Y,rot90(DT_prom,3),10), colormap('jet')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'DT [-]';
c.Label.FontSize = 16;	
c.LineWidth = 1.2;
c.Limits = [round(min(min(DT_prom)),1) round(max(max(DT_prom)),1)];
filename2 = strcat(filename,'-',caso,caso2);
saveas(fig,fullfile(fpath,filename2),'epsc');