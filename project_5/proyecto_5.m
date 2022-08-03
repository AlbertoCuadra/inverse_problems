%% Proyecto 5 - Reconstrucción Algebraica
clear all; clc; close all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
clear red blue green purple yellow orange;
% ------------------------------------
tomografia = '1';
load(strcat('DatosTomografia',tomografia,'.mat'));
% A == imagen 70*70 pixels
% R == sinograma correspondiente a 50 angulos
% entre 1 y 180 grados
F= A\R(:);
aux = strcat('DatosTomografia',tomografia,'.mat');
save(strcat('Fdatos',tomografia));
%% Tomografia
clear all; clc; close all;
tomografia = '4';
load(strcat('Fdatos',tomografia));
n = length(F(:,1));
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
grid on; box on;hold on; axis tight
pcolor(reshape(F,70,70)); 
% c = colorbar;
% c.LineWidth = 1.2;
% c.Limits = [round(min(min(F)),1) round(max(max(F)),1)];
set(gca,'xTick',[ ],'yTick',[ ])
shading interp;
% colormap ('jet');
colormap ('gray');
saveas(fig,strcat('Figure',tomografia,'1.png'));
%% Sinograma
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',20,'BoxStyle','full','Layer','top')
set(fig,'innerposition',[150,100,900,800],...
    'outerposition',[150 100 900 700]);
box on;hold on; axis tight;grid on;
pcolor(R); 
shading interp;
colormap ('jet');
% colormap ('gray');
c = colorbar;
c.LineWidth = 1.2;
c.Limits = [round(min(min(R)),1) round(max(max(R)),1)];
xlabel('$\theta$','FontSize',24,'interpreter','latex');
ylabel('$t$','FontSize',24,'interpreter','latex');
saveas(fig,strcat('Figure',tomografia,'2.png'));
%% VIDEO Regularizacion por Tikhonov
alpha = [1e-5,1e-4,1e-3,1e-2,0.1,1,10,100,1000,10000];
filename = strcat('tik',tomografia);
fpath = strcat('D:\Google Drive\M2i\Problemas inversos\Proyecto_5\Latex\Figures\tik',...
    tomografia);
mov=VideoWriter(strcat(filename),'MPEG-4');
set(mov,'FrameRate',2);
open(mov);
set(gca,'nextplot','replacechildren');
for i=1:length(alpha)
    selec = alpha(i);
    F_new(1,:) = (A'*A + selec.*eye(n))\(A'*R(:)); % k=0
    fig = figure;
    set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
    grid on; box on;hold on; axis tight
    pcolor(reshape(F_new,70,70)); 
    set(gca,'xTick',[ ],'yTick',[ ])
    shading interp; 
    colormap ('jet');
    tit = ['$\alpha = $',num2str(selec)];
    title(tit,'interpreter','latex','FontSize',20);
    set(fig,'Visible','off');
    MM(i)=getframe(fig);
    filename2 = strcat(filename,'-',num2str(i),'.png');
    saveas(fig,fullfile(fpath,filename2));
end
writeVideo(mov,MM);
% movie(M);
close(mov);
%% SELECCION DEL ALPHA OPTIMO
selec = 1;
save(strcat('selec',tomografia),'selec');
%% TSVD 
% Calculo de descomposicion en valores singulares (svd)
[U,S,V] = svd(A);% U = M x M, V = N x N
[M,N] = size(S);
save(strcat('SVDdatos',tomografia),'U','S','V','M','N');
%%
clear all; clc; close all;
tomografia = '4';
load(strcat('Fdatos',tomografia));
load(strcat('SVDdatos',tomografia));
r = rank(A);
% Representacion del espectro de los valores singulares
for i=1:r
    x(i) = i;
    sigma(1,i) = S(i,i);
end
% GRAFICA valores singulares
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',14,'BoxStyle','full')
set(fig,'innerposition',[150,100,1200,800],...
    'outerposition',[150 100 1100 700]);
grid on; box on; axis tight; hold on
plot(x(1:1:end),sigma(1:1:end),'d-','color',color(2,:),'LineWidth',1.5,...
  'MarkerFaceColor',color(2,:),'MarkerEdgeColor','black')
set(gca,'yscale','log')
xlabel('$i$','FontSize',20,'interpreter','latex');
ylabel('$\sigma$','FontSize',20,'interpreter','latex');
if max(x)>4000
    xlim([1,4500]);
elseif max(x)>2000
    xlim([1,2500]);
end
saveas(fig,strcat('Figure',tomografia,'4.png'));
%% Video modos
Nmodosf = 2000; 
Deltamodos = 50;
filename = strcat('modos',tomografia);
fpath = strcat('D:\Google Drive\M2i\Problemas inversos\Proyecto_5\Latex\Figures\prueba',...
    tomografia);
mov=VideoWriter(strcat(filename),'MPEG-4');
set(mov,'FrameRate',10);
open(mov);
set(gca,'nextplot','replacechildren');
for i=1:round(Nmodosf/Deltamodos)
    m = i*Deltamodos; %Modos!
%     B = U(:,1:m)*S(1:m,1:m)*V(:,1:m)';
%     F_svd= B'*R(:);
%     B = V(:,1:m)*inv(S(1:m,1:m))*U(:,1:m)';
    B = V(:,1:m)*pinv(S(1:m,1:m))*U(:,1:m)';
    F_svd= B*R(:);
    fig = figure;
    set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
    grid on; box on;hold on; axis tight;
    pcolor(reshape(F_svd,70,70)); 
    shading interp;
    set(gca,'xTick',[ ],'yTick',[ ])
    colormap ('jet');
    tit = ['modos $=\ $',num2str(x(i*Deltamodos))];
    title(tit,'interpreter','latex','FontSize',20);
    set(fig,'Visible','off');
    MM(i)=getframe(fig);
    filename2 = strcat(filename,'-',num2str(i),'.png');
    saveas(fig,fullfile(fpath,filename2));
end
writeVideo(mov,MM);
% movie(M);
close(mov);
%% DEFECTOS
tomografia='3';
%% DEFECTOS Tikhonov
imagen=strcat('tik',tomografia,'-5','.png');
filename = strcat('tik',tomografia);
fpath = strcat('D:\Google Drive\M2i\Problemas inversos\Proyecto_5\Latex\Figures\tik',...
    tomografia,'\',imagen);
imagen = importdata(fpath);
%  Threshold
[imagenBW,imagenMasked] = createMask2(imagen);
fig = imshow(imagenMasked,'border','tight');
saveas(fig,strcat('Figure',tomografia,'5.png'));
%% DEFECTOS TSVD
imagen=strcat('modos',tomografia,'-20','.png');
filename = strcat('modo',tomografia);
fpath = strcat('D:\Google Drive\M2i\Problemas inversos\Proyecto_5\Latex\Figures\modo',...
    tomografia,'\',imagen);
imagen = importdata(fpath);
% Threshold
[imagenBW,imagenMasked] = createMask2(imagen);
fig = imshow(imagenMasked,'border','tight');
saveas(fig,strcat('Figure',tomografia,'6.png'));