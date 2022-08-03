%% Proyecto 5 - Reconstrucción Algebraica
clear all; clc; close all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
clear red blue green purple yellow orange;
% ------------------------------------
load('DatosTomografia2.mat');
% A == imagen 70*70 pixels
% R == sinograma correspondiente a 50 angulos
% entre 1 y 180 grados
F= A\R(:);
save('Fdatos2');
%% Tomografia
clear all; clc; close all;
load('Fdatos2');
n = length(F(:,1));
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
grid on; box on;hold on; axis tight
pcolor(reshape(F,70,70)); 
% c = colorbar;
% c.LineWidth = 1.2;
% c.Limits = [round(min(min(F)),1) round(max(max(F)),1)];
set(gca,'xTick',[ ],'yTick',[ ])
shading interp;
% colormap ('jet');
colormap ('gray');
saveas(fig,'Figure21.eps');
%% Sinograma
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
set(fig,'innerposition',[150,100,900,800],...
    'outerposition',[150 100 900 700]);
box on;hold on; axis tight;grid on;
pcolor(R); 
shading interp;
colormap ('jet');
% colormap ('gray');
c = colorbar;
c.LineWidth = 1.2;
c.Limits = [round(min(min(R)),1) round(max(max(R)),1)];
xlabel('$\theta$','FontSize',16,'interpreter','latex');
ylabel('$t$','FontSize',16,'interpreter','latex');
saveas(fig,'Figure22.eps');
%% Regularizacion por Tikhonov
alpha = [1e-5,1e-3,0.1,1,10,100,1000];
selec = alpha(6);
F_new(1,:) = (A'*A + selec.*eye(n))\(A'*R(:)); % k=0
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
grid on; box on;hold on; axis tight
pcolor(reshape(F_new,70,70)'); 
% c = colorbar;
% c.LineWidth = 1.2;
% c.Limits = [round(min(min(F_new)),1) round(max(max(F_new)),1)];
set(gca,'xTick',[ ],'yTick',[ ])
shading interp; 
colormap ('jet');
tit = ['$\alpha = $',num2str(selec)];
title(tit,'interpreter','latex','FontSize',20);
saveas(fig,'Figure23.eps');
%% SELECCION DEL ALPHA OPTIMO
selec = 1e-1;
save('selec1','selec');
%% TSVD 
% Calculo de descomposicion en valores singulares (svd)
[U,S,V] = svd(A);% U = M x M, V = N x N
[M,N] = size(S);
% save('SVDdatos1','U','S','V','M','N');
%%
clear all; clc; close all;
load('Fdatos1');
load('SVDdatos1');
r = rank(A);
% Representacion del espectro de los valores singulares
for i=1:r
    x(i) = i;
    sigma(1,i) = S(i,i);
end
% GRAFICA valores singulares
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1200,800],...
    'outerposition',[150 100 1100 700]);
grid on; box on; axis tight; hold on
plot(x(1:30:end),sigma(1:30:end),'d-','color',color(2,:),'LineWidth',1.5,...
  'MarkerFaceColor',color(2,:),'MarkerEdgeColor','black')
set(gca,'yscale','log')
xlabel('$i$','FontSize',20,'interpreter','latex');
ylabel('$\sigma$','FontSize',20,'interpreter','latex');
xlim([1,4500]);
saveas(fig,'Figure24.eps');
%% 
m = 100; %Modos!
B = U(:,1:m)*S(1:m,1:m)*V(:,1:m)';
F_svd= B'*R(:);
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
grid on; box on;hold on; axis tight;
pcolor(reshape(F_svd,70,70)); 
shading interp;
% c = colorbar;
% c.LineWidth = 1.2;
% c.Limits = [round(min(min(F_svd)),1) round(max(max(F_svd)),1)];
set(gca,'xTick',[ ],'yTick',[ ])
% colormap ('gray');
colormap ('jet');
saveas(fig,'Figure25.eps');