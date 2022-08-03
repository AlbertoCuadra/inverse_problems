%% Proyecto 2 - Crimenes inversos
%% Parte 3 - Diferencias regresivas d = -d
clear all; close all; clc;
% load('datos_proyecto2');
load('datos_proyecto22');
k = 9;
A = full(gallery('tridiag',length(x),-d,1+2*d,-d)); % Parte 3
% RUIDO
eps = 0;
u_dif(:,end) = u_dif(:,end) + eps*rand(size(u_dif(:,end)));
% u_dif(:,end) = u_dif(:,end) + eps*ones(size(u_dif(:,end)));
% Para el apartado 2.3, partimos de la solucion exacta con ruido
t = 0.4;
u_exacta = 10*exp(-4*t)*sin(2*x(1:end));
u_exacta = u_exacta';
u_dif2(:,1) = u_exacta(:,end) + eps*rand(size(u_exacta(:,end)));
%  u_dif2(:,1) = u_exacta(:,end) + eps*ones(size(u_exacta(:,end)));
% -------------------------------------------------------
% Calculamos u en un instante T_0 a partir de mediciones en un instante
% posterior T
u_new(:,1) = u_dif(:,end);% Para tener como primer vector la solucion 
% en el instante T_final
for i=1:k
    u_new(:,i+1) = A*u_new(:,i); % Apartado 3
end
t = 0.4-k*deltat;
u_exacta = 10*exp(-4*t)*sin(2*x(1:end));
u_exacta = u_exacta';

% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
plot(x,u_exacta,'color',color(2,:),'LineWidth',1.5)
plot(x,u_new(:,end),'d','color',color(3,:),'LineWidth',1.5,...
    'MarkerFaceColor',color(3,:),'MarkerEdgeColor','black')
% REPRESENTACION DE RESULTADOS
xlabel('$x$','FontSize',20,'interpreter','latex');
ylabel('$u(x,t)$','FontSize',20,'interpreter','latex');
tit = ['$\Delta x = $',num2str(round(deltax,4)),...
    ', $\Delta t = $',num2str(round(deltat,4))];
tit2 = ['$T = $',num2str(round(t,3)),...
    ', $\epsilon = $',num2str(round(eps,4))];
title({tit,tit2},'interpreter','latex');
leg{1} = ['Exacta'];
leg{2} = ['Aproximaci\''on'];
legend(leg,'FontSize',16,'Location','northeast','interpreter','latex')

error = norm(u_exacta(:,1)-u_new(:,end))/norm(u_exacta(:,end));
error_ab = norm(u_exacta(:,1)-u_new(:,end));
fprintf('Error relativo = %.3e \n',error)
fprintf('Error absoluto = %.3e \n',error_ab)
% Metodo muy sensible a ruido
% Crimenes inversos -> utilizar el mismo metodo la inversa da lugar a
% mejores resultados, pero es falso, a esto se le conoce como crimenes
% inversos.