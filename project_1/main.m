% PROYECTO 1 - PROSPECCION GRAVITATORIA
% Realizado por Alberto Cuadra Lara
% Ultima mod. 16/02/2018
help proyecto_1
clc; clear all; close all;
global color
% ---------------------------------------------------------
% Colores
red = [0.85,0.33,0.1]; blue=[0,0.45,0.74]; green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
% ---------------------------------------------------------
% CONFIGURACION DEL PROBLEMA
% Lectura de datos
load('CampoGravitacional2.mat');
t = F(1,:);
eps = 0.0001;
B = F(2,:) + eps*rand(size(t)); %vector independiente.
% GRAFICA DATOS F_z (equivalente a B)
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[250,200,1150,700],...
    'outerposition',[250 200 1000 600]);
hold on; grid on; grid minor; box on; axis tight;
plot(t,B,'d-','color',color(2,:),'LineWidth',1.5,...
    'MarkerFaceColor',color(3,:),'MarkerEdgeColor','black'); 
xlabel('$t$','FontSize',20,'interpreter','latex');
ylabel('$F_{z,2}$','FontSize',20,'interpreter','latex');
leg{1} = ['$F_{z,2}$'];
% tit = ['$\epsilon = $',num2str(eps)];
% title(tit,'interpreter','latex');
legend(leg,'FontSize',16,'Location','northeast','interpreter','latex')
% ---------------------------------------------------------
G = -1; % Simplificacion constante gravitacional universal
deltax = 0.01; 
x = 0:deltax:1; % Discretizacion espacial
% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[250,200,1150,700],...
    'outerposition',[250 200 1000 600]);
% ---------------------------------------------------------
for k=1:5 % Grado del polinomio
    for i=1:length(t)
        f = ((x-t(i)).^2+1).^(-3/2);
        A(i,1) = -G*trapz(x,f);
        for j=2:k+1
            f = x.^(j-1)./((x-t(i)).^2+1).^(3/2);
            A(i,j) = -G*trapz(x,f);
        end
    end
    % Resolvemos el problema de minimizacion por minimos cuadrados -> X
    X = (transpose(A)*A)\(transpose(A)*B'); % coeficientes a_i
    rho(k) = cond(A);
    % GRAFICA SOLUCION
    hold on; grid on; box on; axis tight;
    plot(x,polyval(flip(X),x),'LineWidth',1.5,'color',color(k,:))
    xlabel('$x$','FontSize',20,'interpreter','latex');
    ylabel('$\lambda$','FontSize',20,'interpreter','latex');
    leg{k} = ['Grado ',num2str(k)];
    tit = ['$\epsilon = $',num2str(eps)];
    title(tit,'interpreter','latex');
    legend(leg,'FontSize',16,'Location','northwest','interpreter','latex')
end
%% Comprobacion del polinomio
lambda = polyval(flip(X),x);
% Configuracion grafica
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[250,200,1150,700],...
    'outerposition',[250 200 1000 600]);
% ---------------------------------------------------------
for i=1:length(t)
    f = lambda./((x-t(i)).^2+1).^(3/2);
    F_z(i) = -G*trapz(x,f);
    error(i) = (norm(B(i)-F_z(i)))/(1+norm(F_z(i)))*100;
end
hold on; grid on; box on; axis tight;
plot(t,error,'LineWidth',1.5,'color',color(2,:))
xlabel('$t$','FontSize',20,'interpreter','latex');
ylabel('$\epsilon\ [\%]$','FontSize',20,'interpreter','latex');
