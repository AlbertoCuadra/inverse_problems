%% Proyecto 2 - Crimenes inversos
%% Parte 1
clc; clear all; close all
% Colores
red = [0.85,0.33,0.1]; blue=[0,0.45,0.74]; green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56]; yellow = [0.93,0.69,0.13];
% ---------------------------------------------------------------
% PARAMETROS DEL PROBLEMA
L = pi;% Longitud de la barra
T = 0.4;% Tiempo final
M = 40;
N = 400;
deltax = L/(M+1);
deltat = T/(N+1);
D = 1; % Coef. de difusion
for i=1:M
    x(i) = i*deltax;% Discretizacion espacial
end
for i=1:N
    t(i) = i*deltat;% Discretizacion temporal
end
% Condiciones iniciales
f = 10.*sin(2*x);
% DIFERENCIAS FINITAS u(x,t)
d = D*deltat/(deltax)^2;
A = full(gallery('tridiag',length(x),d,1-2*d,d));
u(:,1) = f;
for i=1:length(t)+1
    u(:,i+1) = A*u(:,i);
end

% surf(t,x,u)
% Sol. exacta
% x=[0,x,pi];
 t= [0,t,T];
for i=1:length(x)
    for j=1:length(t)
        u_exacta(i,j) = 10*exp(-4*t(j))*sin(2*x(i));
    end
end
% figure
% % surf(t,x,u_exacta)
% ERROR
% error = max(abs(u_exacta(:,end)-u(2:end-1,end)))./max(abs(u_exacta(:,end)));
for i=2:M
    errorn(i) = abs(u_exacta(i,N+2)-u(i,N+2));
    errord(i) = abs(u_exacta(i,N+2));
end
error = max(errorn)/max(errord);

%% Parte 2 
% A es la misma
% RUIDO
eps = 1e-4;
u(:,end) = u(:,end) + eps*rand(size(u(:,end)));
for i=length(t):-1:length(t)-10
    u_new(:,length(t)-i+1) = A\u(:,i);
    u(:,i) = u_new(:,i);
end
t_2 = 0.4-10*deltat;
u_exacta2 = 10*exp(-4*t_2)*sin(2*x(1:end));
u_exacta2 = u_exacta2';
figure
hold on
plot(x,u_new(:,end));
plot(x,u_exacta2);
% Metodo muy sensible a ruido
% Crimenes inversos -> utilizar el mismo metodo la inversa da lugar a
% mejores resultados, pero es falso, a esto se le conoce como crimenes
% inversos.