%% Proyecto 2 - Crimenes inversos
%% Parte 4 - Distinta malla
clear all; close all; clc;
load('datos_proyecto2');
k = 15;
M = 35;
N = 350;
T = 0.4;
for i=1:M% Omitimos extremos (Cond. contorno - conocido)
    x_new(i) = i*L/(M+1);% Discretizacion espacial
end
for i=0:N+1
    t_new(i+1) = i*T/(N+1);% Discretizacion temporal
end
deltat = t_new(2);
u_malla = spline(x,u_dif(:,end),x_new);
A = full(gallery('tridiag',length(x_new),-d,1+2*d,-d)); % Parte 3
% RUIDO
eps = 0;
% u_malla(:,end) = u_malla(:,end) + eps*rand(size(u_malla(:,end)));
u_malla(:,end) = u_malla(:,end) + eps*ones(size(u_malla(:,end)));
u_malla = u_malla';
% Para el apartado 2.3, partimos de la solucion exacta con ruido
t = 0.4;
u_exacta = 10*exp(-4*t)*sin(2*x_new(1:end));
u_exacta = u_exacta';
u_dif2(:,1) = u_exacta(:,end) + eps*rand(size(u_exacta(:,end)));
% -------------------------------------------------------
% Calculamos u en un instante T_0 a partir de mediciones en un instante
% posterior T
u_new = u_malla;% Para tener como primer vector la solucion 
% u_new(:,1) = u_dif2(:,end);% Para tener como primer vector la solucion 
% en el instante T_final
for i=1:k
    u_new(:,i+1) = A*u_new(:,i); % Apartado 3
end
t = 0.4-k*deltat;
u_exacta = 10*exp(-4*t)*sin(2*x(1:end));
u_exacta= spline(x,u_exacta,x_new);
u_exacta = u_exacta';

% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
plot(x_new,u_exacta,'color',color(2,:),'LineWidth',1.5)
plot(x_new,u_new(:,end),'d-','color',color(3,:),'LineWidth',1.5,...
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
% Metodo muy sensible a ruido
% Crimenes inversos -> utilizar el mismo metodo la inversa da lugar a
% mejores resultados, pero es falso, a esto se le conoce como crimenes
% inversos.

error = norm(u_exacta(:,1)-u_new(:,end))/norm(u_exacta(:,end));
error_ab = norm(u_exacta(:,1)-u_new(:,end));
fprintf('Error relativo = %.3e \n',error)
fprintf('Error absoluto = %.3e \n',error_ab)