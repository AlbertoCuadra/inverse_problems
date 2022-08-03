%% Proyecto 4 - Reconstruccion de bases de datos
clear all; clc; close all;
Nn = 20;
x = linspace(0,1,Nn);
y = linspace(0,1,Nn);
z = [0:0.5:1];
puntos = Nn^2;
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
end
%% Formamos matrices con datos incompletos -> A_gap
gaps_porcentaje = 0.6;% En porcentaje [0 1]
gaps = gaps_porcentaje*puntos;
for k=1:length(z)
    i = 1;
    A_gap(:,:,k) = A(:,:,k);
for i=1:gaps
    a = randi(1 size(
    
    
%     nx = randi([2 Nn-1],1);
%     ny = randi([2 Nn-1],1);
%     A_gap(nx,ny,k) = 0;
%     v(i,1,k) = nx; % Almacenamos la posicion de donde variamos los valores
%     v(i,2,k) = ny; % Almacenamos la posicion de donde variamos los valores
%     i = i+1;
end
figure 
surf(y,x,A_gap(:,:,k)); 
ylabel('$y$','FontSize',20,'interpreter','latex');
xlabel('$x$','FontSize',20,'interpreter','latex');
zlabel('$f$','FontSize',20,'interpreter','latex');
end
% %% Reconstruccion de la matriz paso s=0 
% for k=1:length(z)
%     A_recons(:,:,k) = A_gap(:,:,k);
%     %%Rellenamos los huecos haciendo la media
%     for i=1:gaps
%         A_recons(v(i,1,k),v(i,2,k),k) = 0.25*(A_gap(v(i,1,k)+1,v(i,2,k),k)...
%             + A_gap(v(i,1,k)-1,v(i,2,k),k) + A_gap(v(i,1,k),v(i,2,k)+1,k)...
%             + A_gap(v(i,1,k),v(i,2,k)-1,k));
%     end
% figure 
% surf(y,x,A_recons(:,:,k)); 
% end
%% Reconstruccion SVD paso s>=1
k = 1;
maxit = 1000;
A_recons2(:,:,k) = A_gap(:,:,k);
for i=1:maxit
% Calculo de descomposicion en valores singulares (svd)
[U,S,V] = svd(A_recons2(:,:,k));% U = M x M, V = N x N
[M,N] = size(S);
% Representacion del espectro de los valores singulares
    for i=1:N
        x(i) = i;
        sigma(1,i) = S(i,i);
    end
    r = rank(A_gap(:,:,k));
    m = 3; %Modos!
    B = U(:,1:m)*S(1:m,1:m)*V(:,1:m)';
    for i=1:gaps
        A_recons2(v(i,1,k),v(i,2,k),k) = B(v(i,1,k),v(i,2,k));
    end
end
figure 
surf(y,x,A_recons2(:,:,k)); 
RMSE = norm(A(:,:,k)-A_recons2(:,:,k))/Nn
E = max(max(abs(A(:,:,k)-A_recons2(:,:,k))))