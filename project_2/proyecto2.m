%% Proyecto 2 - Crimenes inversos
%% Parte 1 - Problema directo
clc; clear all; close all
global color
format shortEng
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
% ---------------------------------------------------------------
% PARAMETROS DEL PROBLEMA
L = pi;% Longitud de la barra
T = 0.385;% Tiempo final
M  = 35;% Numero de nodos de la discretizacion espacial
N = 350; % Numero de intervalos temporales
D = 1; % Coef. de difusion
% DISCRETIZACION ESPACIAL Y TEMPORAL
for i=1:M% Omitimos extremos (Cond. contorno - conocido)
    x(i) = i*L/(M+1);% Discretizacion espacial
end
for i=0:N+1
    t(i+1) = i*T/(N+1);% Discretizacion temporal
end
% CONDICION INICIAL
f = 10.*sin(2*x);
% CALCULO POR DIFERENCIAS FINITAS u(x,t)
deltat = t(2);
deltax = x(1);
d = D*deltat/(deltax)^2;
A = full(gallery('tridiag',length(x),d,1-2*d,d));
u_dif(:,1) = f;
for i=1:N+1
    u_dif(:,i+1) = A*u_dif(:,i);
end
% CALCULO DE LA SOLUCION EXACTA u(x,t)
for i=1:length(x)
    for j=1:length(t)
        u_exacta(i,j) = 10*exp(-4*t(j))*sin(2*x(i));
    end
end
% CALCULO DEL ERROR RELATIVO EN EL INSTANTE FINAL T
error = norm(u_exacta(:,end)-u_dif(:,N+2))/norm(u_exacta(:,end));
% REPRESENTACION DE LOS RESULTADOS
surf(t,x,u_exacta,'LineStyle','none');
xlabel('$t$','FontSize',20,'interpreter','latex');
ylabel('$x$','FontSize',20,'interpreter','latex');
zlabel('$u(x,t)$','FontSize',20,'interpreter','latex');
% plot(x,u_dif(:,end));
% IMPRESION DE RESULTADOS
fprintf('T = %.2f \n',T)
fprintf('M = %d \n',M)
fprintf('N = %d \n',N)
fprintf('Deltax = %.4f \n',deltax)
fprintf('Deltat = %.3e \n',deltat)
fprintf('d = D*Deltat/Deltax^2 = %.4f \n',d)
fprintf('Error = %.3e \n',error)
%% Parte 2 - Problema inverso
clear all; close all; clc;
% load('datos_proyecto2');
load('datos_proyecto22');
k = 10;
% RUIDO
eps = 0;
% u_dif(:,end) = u_dif(:,end) + eps*rand(size(u_dif(:,end)));
u_dif(:,end) = u_dif(:,end) + eps*ones(size(u_dif(:,end)));
% Para el apartado 2.3, partimos de la solucion exacta con ruido
t = 0.4;
u_exacta = 10*exp(-4*t)*sin(2*x(1:end));
u_exacta = u_exacta';
u_dif2(:,1) = u_exacta(:,end) + eps*rand(size(u_exacta(:,end)));
%
% -------------------------------------------------------
% Calculamos u en un instante T_0 a partir de mediciones en un instante
% posterior T
u_new(:,1) = u_dif(:,end);% Para tener como primer vector la solucion 
% en el instante T_final
for i=1:k
    u_new(:,i+1) = A\u_new(:,i);
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
plot(x,u_new(:,end),'d-','color',color(3,:),'LineWidth',1.5,...
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