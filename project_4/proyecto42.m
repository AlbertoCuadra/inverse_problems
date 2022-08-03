%% Proyecto 4 - Reconstruccion de bases de datos
clear all; clc; close all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
% ---------------------------------------------------------------- %
N = 50;
x = linspace(0,1,N);
y = linspace(0,1,N);
z = [0:0.5:1];
puntos = N^2;
for i=1:length(z)
    for j=1:length(x)
        for k=1:length(y)
            aux(j,k,i) = f(x(j),y(k),z(i));
        end
    end
F = 2.*(aux(:,:,i)-min(min(aux(:,:,i))))./(max(max(aux(:,:,i)))-min(min(aux(:,:,i))))-1; 
% Funcion normalizada con valores entre 0 y 1
% Para un valor fijo de z = [0,1], construimos una base de datos 
A(:,:,i) = F;
figure
surf(y,x,A(:,:,i)); 
ylabel('$y$','FontSize',20,'interpreter','latex');
xlabel('$x$','FontSize',20,'interpreter','latex');
zlabel('$f$','FontSize',20,'interpreter','latex');
    if i==1
        tit = ['Matriz original $A^1$',': ',...
            '$z = $',num2str(z(i)),', ',...
            '$N = $',num2str(N)];
        title(tit,'interpreter','latex','FontSize',16);
    elseif i==2
        tit = ['Matriz original $A^2$',': ',...
            '$z = $',num2str(z(i)),', ',...
            '$N = $',num2str(N)];
        title(tit,'interpreter','latex','FontSize',16);
    else
        tit = ['Matriz original $A^3$',': ',...
            '$z = $',num2str(z(i)),', ',...
            '$N = $',num2str(N)];
        title(tit,'interpreter','latex','FontSize',16);
    end
end
%% Formamos matrices con datos incompletos -> A_gap
gaps_porcentaje = 0.6;% En porcentaje [0 1]
gaps = gaps_porcentaje*puntos;
for k=1:length(z)
    [X,Y] = meshgrid(1:N,1:N);
    x_r = X(:);
    y_r = Y(:);
    x_r = [x_r,y_r];
    A_gap(:,:,k) = A(:,:,k);
for i=1:gaps
    a = randi([1 size(x_r,1)],1);
    nx = x_r(a,1);
    ny = x_r(a,2);
    x_r(a,:) = [];
%     A_gap(nx,ny,k)=0;
    if (nx==N && ny ==N)
        A_gap(nx,ny,k)=(A(nx-1,ny,k) + A(nx,ny-1,k))/2;
    elseif (nx==1 && ny==1)
        A_gap(nx,ny,k)=(A(nx+1,ny,k) + A(nx,ny+1,k))/2;
    elseif (nx==1 && ny==N)
        A_gap(nx,ny,k)=(A(nx+1,ny,k) + A(nx,ny-1,k))/2;
    elseif (nx==N && ny==1)
        A_gap(nx,ny,k)=(A(nx-1,ny,k) + A(nx,ny+1,k))/2;
    elseif (nx<N && ny<N && nx>1 && ny>1)
        A_gap(nx,ny,k)=(A(nx+1,ny,k) + A(nx-1,ny,k) + A(nx,ny+1,k) + A(nx,ny-1,k))/4;
    elseif (nx==N && ny <N)
        A_gap(nx,ny,k)=(A(nx-1,ny,k) + A(nx,ny+1,k) + A(nx,ny-1,k))/3;
    elseif (nx==1 && ny <N)
        A_gap(nx,ny,k)=(A(nx+1,ny,k) + A(nx,ny+1,k) + A(nx,ny-1,k))/3;
    elseif (nx<N && ny==N)
        A_gap(nx,ny,k)=(A(nx+1,ny,k) + A(nx-1,ny,k) + A(nx,ny-1,k))/3;
    elseif (nx<N && ny==1)
        A_gap(nx,ny,k)=(A(nx+1,ny,k) + A(nx-1,ny,k) + A(nx,ny+1,k))/3;
    end
    v(i,1,k) = nx; % Almacenamos la posicion de donde variamos los valores
    v(i,2,k) = ny; % Almacenamos la posicion de donde variamos los valores
end
figure 
surf(y,x,A_gap(:,:,k)); 
ylabel('$y$','FontSize',20,'interpreter','latex');
xlabel('$x$','FontSize',20,'interpreter','latex');
zlabel('$f$','FontSize',20,'interpreter','latex');
if k==1
    tit = ['Matriz gappy $A^1$',': ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',16);
elseif k==2
    tit = ['Matriz gappy $A^2$',': ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',16);
else
    tit = ['Matriz gappy $A^3$',': ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',16);
    end
end
%% Reconstruccion SVD paso s>=1
k = 1;
maxit = 2000;
A_recons2(:,:,k) = A_gap(:,:,k);

if k==1
    aux=0;
elseif k==2
    aux=0.5;
else
    aux=1;
end

for m=1:6 %Modos!
for i=1:maxit
% Calculo de descomposicion en valores singulares (svd)
[U,S,V] = svd(A_recons2(:,:,k));% U = M x M, V = N x N
[M,N] = size(S);
% Representacion del espectro de los valores singulares
    for i=1:N
        x_espectro(i) = i;
        sigma(1,i) = S(i,i);
    end
    r = rank(A_gap(:,:,k));
%     m = 3; 
    B = U(:,1:m)*S(1:m,1:m)*V(:,1:m)';
    for i=1:gaps
        A_recons2(v(i,1,k),v(i,2,k),k) = B(v(i,1,k),v(i,2,k));
    end
end
figure 
axis tight;
surf(y,x,A_recons2(:,:,k)); 
xlabel('$x$','FontSize',20,'interpreter','latex');
ylabel('$y$','FontSize',20,'interpreter','latex');
zlabel('$f$','FontSize',20,'interpreter','latex');
if k==1
    tit = ['Matriz reconstruida $A^1$',': ',...
        '$m = $',num2str(m),', ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',16);
elseif k==2
    tit = ['Matriz reconstruida $A^2$',': ',...
        '$m = $',num2str(m),', ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',16);
else
    tit = ['Matriz reconstruida $A^3$',': ',...
        '$m = $',num2str(m),', ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',16);
end

RMSE(m) = norm(A(:,:,k)-A_recons2(:,:,k))/N;
E(m) = max(max(abs(A(:,:,k)-A_recons2(:,:,k))));
% fprintf('RMSE = %.3e \n',RMSE(m))
% fprintf('MaxE = %.3e \n',E(m))
end
%% REPRESENTACION DEL ERROR
m = 1:1:20;
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',18,'BoxStyle','full')
set(fig,'innerposition',[150,100,900,800],...
    'outerposition',[150 100 900 700]);
grid on; box on;hold on; axis tight
plot(m,RMSE,'LineWidth',2,'color',color(2,:));
plot(m,E,'LineWidth',2,'color',color(3,:));
set(gca,'yscale','log')
xlabel('$m$','FontSize',20,'interpreter','latex');
ylabel('E','FontSize',20,'interpreter','latex');
leg{1} = ['RMSE'];
leg{2} = ['maxE'];
legend(leg,'FontSize',22,'Location','northeast','interpreter','latex')
if k==1
    tit = ['Error matriz $A^1$: $N = $',num2str(N),', ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',26);
elseif k==2
    tit = ['Error matriz $A^2$: $N = $',num2str(N),', ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',26);
else
    tit = ['Error matriz $A^3$: $N = $',num2str(N),', ',...
        num2str(gaps_porcentaje*100),'$\%$ gappyness'];
    title(tit,'interpreter','latex','FontSize',26);
end