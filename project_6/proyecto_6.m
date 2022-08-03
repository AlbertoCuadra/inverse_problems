%% Proyecto 6 - Derivadas topologicas
clear all; clc; close all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
clear red blue green purple yellow orange
% ------------------------------------
load('DatosDT.mat');
N = 50;
x_mesh = linspace(-1,1,N);
y_mesh = linspace(-1,1,N);
[X,Y] = meshgrid(x_mesh,y_mesh);
x = X(:);
y = Y(:);
x = [x,y];
DT = zeros(length(Angulos),length(x));
w = zeros(length(Angulos),length(x));
DT_new = zeros(1,length(x));
OndaMedida = OndaMedida.';% .' para solo transponer
[Nr,Nd] = size(OndaMedida);
for j=1:length(Xreceptores)
    u_incidente_receptores(:,j) = exp(i*lambda.*Angulos.*Xreceptores(j).*Yreceptores(j));
end
d = [cos(Angulos)',sin(Angulos)'];
for k=1:length(Angulos)
    u_incidente(k,:) = exp(i*lambda*x*d(k,:)');
end
for k=1:length(Angulos)
    for j=1:length(x)
        for t=1:length(Xreceptores)
            z = lambda*norm(x(j)-Xreceptores(t).*Yreceptores(t));
            w(k,t) = w(k,t) + (i/4)*besselh(0,1,z)*conj(OndaMedida(k,t)-u_incidente_receptores(k,t));
        end
        DT(k,j) = real(u_incidente(k,t).*w(k,t));
    end
    DT_new = DT_new + DT(k,:);
end
DT_new = reshape(DT_new,[N,N]);
pcolor(x_mesh,y_mesh,DT_new)