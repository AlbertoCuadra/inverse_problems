%% Proyecto 2 - Crimenes inversos
%% Parte 1 - Convergencia Variando N
clc; clear all; close all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
% ---------------------------------------------------------------
% PARAMETROS DEL PROBLEMA
L = pi;% Longitud de la barra
T = 0.1;% Tiempo final
M  = 20;% Numero de nodos de la discretizacion espacial
Nf = 50; % Numero de intervalos temporales
D = 1; % Coef. de difusion
% DISCRETIZACION ESPACIAL Y TEMPORAL

lambda = 1;
% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
for T=2:0.2:2
    k = 1;
for N=20:Nf
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
error(k) = norm(u_exacta(:,end)-u_dif(:,N+2))/norm(u_exacta(:,end));
% REPRESENTACION DE LOS RESULTADOS
% surf(t,x,u_dif);
% plot(x,u_dif(:,end));
% IMPRESION DE RESULTADOS
% fprintf('T = %.2f \n',T)
% fprintf('M = %d \n',M)
% fprintf('N = %d \n',N)
% fprintf('Deltax = %.4f \n',deltax)
% fprintf('Deltat = %.3e \n',deltat)
% fprintf('d = D*Deltat/Deltax^2 = %.4f \n',d)
% fprintf('Error = %.3e \n',error)
d_vector(k) = d;
k = k+1;
clear u_dif u_exacta x t
end
N=[20:Nf];
plot(N,error,'LineWidth',1.5)%Error en funcion de variar N
% plot(N,d_vector,'LineWidth',1.5)
% REPRESENTACION DE RESULTADOS
xlabel('$N$','FontSize',20,'interpreter','latex');
% ylabel('$d$','FontSize',20,'interpreter','latex');
ylabel('$\epsilon$','FontSize',20,'interpreter','latex');
leg{lambda} = ['$d(N)$ con $M=20$ y $T = $',num2str(T)];
lambda = lambda +1;
clear error
end
% plot(N,ones(length(N))*0.5,'color',color(3,:),'LineWidth',1.5)
% leg{lambda} = ['L\''imite convergencia'];
legend(leg,'FontSize',16,'Location','northeast','interpreter','latex')
%% Parte 1 - Convergencia Variando M
clc; clear all; close all
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
% ---------------------------------------------------------------
% PARAMETROS DEL PROBLEMA
L = pi;% Longitud de la barra
T = 0.4;% Tiempo final
Mf  = 40;% Numero de nodos de la discretizacion espacial
N = 500; % Numero de intervalos temporales
D = 1; % Coef. de difusion
lambda = 1;
% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
for T=0.4:0.2:1
    k = 1;
for M=20:Mf
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
error(k) = norm(u_exacta(:,end)-u_dif(:,N+2))/norm(u_exacta(:,end));
% REPRESENTACION DE LOS RESULTADOS
% surf(t,x,u_dif);
% plot(x,u_dif(:,end));
% IMPRESION DE RESULTADOS
% fprintf('T = %.2f \n',T)
% fprintf('M = %d \n',M)
% fprintf('N = %d \n',N)
% fprintf('Deltax = %.4f \n',deltax)
% fprintf('Deltat = %.3e \n',deltat)
% fprintf('d = D*Deltat/Deltax^2 = %.4f \n',d)
% fprintf('Error = %.3e \n',error)
d_vector(k) = d;
k = k+1;
clear u_dif u_exacta x t
end
M=[20:Mf];
plot(M,error,'LineWidth',1.5)%Error en funcion de variar N
% plot(M,d_vector,'LineWidth',1.5)
% REPRESENTACION DE RESULTADOS
xlabel('$M$','FontSize',20,'interpreter','latex');
% ylabel('$d$','FontSize',20,'interpreter','latex');
ylabel('$\epsilon$','FontSize',20,'interpreter','latex');
leg{lambda} = ['$d(M)$ con $N = $',num2str(N),' y $T = $',num2str(T)];
lambda = lambda +1;
clear error
end
% plot(M,ones(length(M))*0.5,'color',color(3,:),'LineWidth',1.5)
% leg{lambda} = ['L\''imite convergencia'];
legend(leg,'FontSize',16,'Location','northeast','interpreter','latex')
%% Parte 1 - Convergencia Variando T
clc; clear all; close all
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
% ---------------------------------------------------------------
% PARAMETROS DEL PROBLEMA
L = pi;% Longitud de la barra
Tf = 1.4;% Tiempo final
M  = 40;% Numero de nodos de la discretizacion espacial
N = 400; % Numero de intervalos temporales
D = 1; % Coef. de difusion
% DISCRETIZACION ESPACIAL Y TEMPORAL
k = 1;
for T=0.4:0.1:Tf
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
error(k) = norm(u_exacta(:,end)-u_dif(:,N+2))/norm(u_exacta(:,end));
% REPRESENTACION DE LOS RESULTADOS
% surf(t,x,u_dif);
% plot(x,u_dif(:,end));
% IMPRESION DE RESULTADOS
% fprintf('T = %.2f \n',T)
% fprintf('M = %d \n',M)
% fprintf('N = %d \n',N)
% fprintf('Deltax = %.4f \n',deltax)
% fprintf('Deltat = %.3e \n',deltat)
% fprintf('d = D*Deltat/Deltax^2 = %.4f \n',d)
% fprintf('Error = %.3e \n',error)
d_vector(k) = d;
k = k+1;
clear u_dif u_exacta x t
end
T=[0.4:0.1:Tf];
% plot(T,error)%Error en funcion de variar N
% CONFIGURACION GRAFICA
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
plot(T,d_vector,'color',color(2,:),'LineWidth',1.5)
plot(T,ones(length(T))*0.5,'color',color(1,:),'LineWidth',1.5)
% REPRESENTACION DE RESULTADOS
xlabel('$T$','FontSize',20,'interpreter','latex');
ylabel('$d$','FontSize',20,'interpreter','latex');
% ylabel('$u(x,t)$','FontSize',20,'interpreter','latex');
leg{1} = ['$d(T)$ con $N=40$ y $N=400$'];
leg{2} = ['L\''imite convergencia'];
legend(leg,'FontSize',16,'Location','northwest','interpreter','latex')
%% Parte 1 - Convergencia Variando T para d fijo
clc; clear all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
% ---------------------------------------------------------------
% PARAMETROS DEL PROBLEMA
L = pi;% Longitud de la barra
Tf = 3.2;% Tiempo final
M  = 40;% Numero de nodos de la discretizacion espacial
D = 1; % Coef. de difusion
% DISCRETIZACION ESPACIAL Y TEMPORAL
df = 0.3;
lambda=1;
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',16,'BoxStyle','full')
set(fig,'innerposition',[150,100,1250,800],...
    'outerposition',[150 100 1100 700]);
hold on; grid on; grid minor; box on; axis tight;
for T=0.4:0.4:Tf
k = 1;
for d=0.05:0.01:df
    
for i=1:M% Omitimos extremos (Cond. contorno - conocido)
    x(i) = i*L/(M+1);% Discretizacion espacial
end
N = round(T/(d*(L/(M+1)^2))-1);
deltax = x(1);
for i=0:N+1
    t(i+1) = i*d*(L/(M+1))^2;% Discretizacion temporal
end
deltat = t(2);
% CONDICION INICIAL
f = 10.*sin(2*x);
% CALCULO POR DIFERENCIAS FINITAS u(x,t)

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
error(k) = norm(u_exacta(:,end)-u_dif(:,N+2))/norm(u_exacta(:,end));
% REPRESENTACION DE LOS RESULTADOS
% surf(t,x,u_dif);
% plot(x,u_dif(:,end));
% IMPRESION DE RESULTADOS
% fprintf('T = %.2f \n',T)
% fprintf('M = %d \n',M)
% fprintf('N = %d \n',N)
% fprintf('Deltax = %.4f \n',deltax)
% fprintf('Deltat = %.3e \n',deltat)
% fprintf('d = D*Deltat/Deltax^2 = %.4f \n',d)
% fprintf('Error = %.3e \n',error)

k = k+1;
clear u_dif t
end
d=0.05:0.01:df;

% CONFIGURACION GRAFICA
plot(d,error,'LineWidth',1.5)
% REPRESENTACION DE RESULTADOS
xlabel('$d$','FontSize',20,'interpreter','latex');
ylabel('$\epsilon$','FontSize',20,'interpreter','latex');
leg{lambda} = ['$T$ = ',num2str(T)];
legend(leg,'FontSize',20,'Location','eastoutside','interpreter','latex')
lambda=lambda +1;
end