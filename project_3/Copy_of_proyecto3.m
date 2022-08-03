%% Proyecto 3 - Ecuacion del calor hacia atras
%% Parte 1 - PROBLEMA DIRECTO
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
T = 0.4;% Tiempo final
M  = 30;% Numero de nodos de la discretizacion espacial
N = 10; %Numero de terminos a considerar en la serie
D = 1; % Coef. de difusion
for i=0:M
    x(i+1) = i*L/M;% Discretizacion espacial
end
delta_x = x(2)-x(1)
y = x; 
f = 2/pi*x.*(x<=pi/2) + 2/pi*(pi-x).*(x>pi/2);
% plot(x,f); %Grafica f (temperatura inicial)
for i=1:length(x)
    for j=1:length(y)
        K(i,j) = funcionK(x(i),y(j),L,T,N,D);
    end
end       
A = L/M*K;
g = A*f';
% plot(x,g) %Grafica g (temperatura final)
% Parte 1- diferencias
% PARAMETROS DEL PROBLEMA
M = 40;
N = 1001*T;
deltax_dif = L/(M+1);
deltat_dif = T/(N+1);
for i=1:M
    x_dif(i) = i*deltax_dif;% Discretizacion espacial
end
for i=1:N
    t_dif(i) = i*deltat_dif;% Discretizacion temporal
end
% Condiciones iniciales
f = 2/pi*x_dif.*(x_dif<=pi/2) + 2/pi*(pi-x_dif).*(x_dif>pi/2);
d = D*deltat_dif/(deltax_dif)^2;
A_dif = full(gallery('tridiag',length(x_dif),d,1-2*d,d));
u_dif(:,1) = f;
for i=1:length(t_dif)+1
    u_dif(:,i+1) = A_dif*u_dif(:,i);
end
uEspectral = spline(x_dif,u_dif(:,end),x);
% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
plot(x,g,'color',color(2,:),'LineWidth',1.5)
plot(x,uEspectral,'d','color',color(3,:),'LineWidth',1.5,...
    'MarkerFaceColor',color(3,:),'MarkerEdgeColor','black')
% REPRESENTACION DE RESULTADOS
xlabel('$x$','FontSize',20,'interpreter','latex');
ylabel('$g(x)$','FontSize',20,'interpreter','latex');
tit = ['$\Delta x = $',num2str(delta_x),', $T = $',num2str(T)];
title(tit,'interpreter','latex');
leg{1} = ['Exacta'];
leg{2} = ['Diferencias'];
legend(leg,'FontSize',16,'Location','northeast','interpreter','latex')
% Calculo del error
error = norm(g' - uEspectral)/norm(g')
%% ALMACENAR DATOS
g = uEspectral;
save('soldirecto');
%% Parte 2 - PROBLEMA INVERSO
clc; clear all; close all
load('soldirecto')
% ---------------------------------------------------------------
% PARAMETROS DEL PROBLEMA
L = pi;% Longitud de la barra
T = 1;% Tiempo final
M  = 30;% Numero de nodos de la discretizacion espacial
N = 10; %Numero de terminos a considerar en la serie
deltax = L/(M+1);
D = 1; % Coef. de difusion
for i=0:M
    x(i+1) = i*L/M;% Discretizacion espacial
end
delta_x = x(2)-x(1);
y = x; 
f = 2/pi*x.*(x<=pi/2) + 2/pi*(pi-x).*(x>pi/2);
for i=1:length(x)
    for j=1:length(y)
        K(i,j) = funcionK(x(i),y(j),L,T,N,D);
    end
end       
A = L/M*K; 
% for i=1:length(x)
%     sum = 0;
%     for j=2:length(x)-1
%         sum = sum + funcionK(x(i),y(j),L,T,N,D);      
%     end
% g(i) = L/(2*M)*(funcionK(x(i),y(1),L,T,N,D) + 2*sum + funcionK(x(i),y(end),L,T,N,D));
% end
% TIKHONOV - CALIBRACION PARAMETRO REGULARIZACION (alpha)
eps = 0;
alpha = [1e-1,4e-2,1e-2,1e-3];
n = length(x);
% Regularizacion de Tikhonov 0,1,2 (k=2)
% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
plot(x,f,'color',color(1,:),'LineWidth',1.5)
leg{1} = ['Exacta'];
% TIKHONOV - CALIBRACION PARAMETRO REGULARIZACION (alpha)
for i=1:length(alpha)
    f_new(i,:) = (A'*A + alpha(i)*eye(n))\(A'*g);
    plot(x,f_new(i,:),'color',color(i+1,:),'LineWidth',1.5)
    leg{i+1} = ['Tikhonov con $\alpha = $',num2str(alpha(i))];
end
% REPRESENTACION DE RESULTADOS
xlabel('$x$','FontSize',20,'interpreter','latex');
ylabel('$f(x)$','FontSize',20,'interpreter','latex');
tit = ['$n = $',num2str(n),', $\epsilon = $',num2str(eps)];
title(tit,'interpreter','latex');
legend(leg,'FontSize',16,'Location','eastoutside','interpreter','latex')
%% TIKHONOV - alpha fijo
eps = 0.001;
g = g + eps*(rand(size(g)));
alpha_f = alpha(2);
f_new(1,:) = (A'*A + alpha_f*eye(n))\(A'*g); % k=0 
% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
plot(x,f,'color',color(1,:),'LineWidth',1.5)
plot(x,f_new(1,:),'color',color(2,:),'LineWidth',1.5)
leg{1} = ['Exacta']; leg{2} = ['Tikhonov 0'];
% METODO ITERATIVO
for i=2:5 % k=1 y k=2, respectivamente
    f_new(i,:) = (A'*A + alpha_f*eye(n))\(A'*g + alpha_f*f_new(i-1,:)');
    plot(x,f_new(i,:),'color',color(i+1,:),'LineWidth',1.5)
    leg{i+1} = ['Tikhonov ',num2str(i-1)];
end
% REPRESENTACION DE RESULTADOS
xlabel('$x$','FontSize',20,'interpreter','latex');
ylabel('$f(x)$','FontSize',20,'interpreter','latex');
tit = ['$n = $',num2str(n),', $\epsilon = $',num2str(eps)];
title(tit,'interpreter','latex');
legend(leg,'FontSize',16,'Location','eastoutside','interpreter','latex')