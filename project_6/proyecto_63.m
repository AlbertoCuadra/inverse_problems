%% Proyecto 6 - Derivadas topológicas
clear all; clc; close all;
global color
% Colores
red = [0.85,0.33,0.1];blue=[0,0.45,0.74];green=[0.47,0.67,0.19];
purple = [0.49,0.18,0.56];yellow = [0.93,0.69,0.13];orange=[0.702,0.349,0];
color=[red;blue;green;purple;yellow;orange];
clear red blue green purple yellow orange
% ------------------------------------
caso='2';
caso2='11';
load(strcat('DatosDTFrec',caso2,'.mat'));
fpath = strcat('D:\Google Drive\M2i\Problemas inversos\Proyecto_6\Latex\Figures');
filename = 'Figure';
stringreceptores = ' 1:20';
N=100;
[X,Y]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));
x=X(:);
y=Y(:);
x=[x,y];
d=[cos(Angulos)' sin(Angulos)'];
[Nr,Nd]=size(OndaMedida);
w=zeros(length(x),Nd);
for k=1:Nd
    u_inc(:,k)=exp(i*lambda*d(k,:)*x');
end
for k=1:Nd  
for a=1:length(x)
    for j=1:1:20
        z=lambda*norm(x(a,:)-[Xreceptores(j) Yreceptores(j)]);
        w(a,k)=w(a,k)+(i/4)*besselh(0,1,z)*conj(OndaMedida(j,k)-exp(lambda*i*[Xreceptores(j) Yreceptores(j)]*d(k,:)'));
    end
    end
end
DT=zeros(N^2,1);
for a=1:length(x)
    for j=1:Nd
        DT(a)=DT(a)+real(u_inc(a,j)*(w(a,j)));
    end
end
MatDT=reshape(DT',N,N);
fig = figure;
set(axes,'LineWidth',1.2,'FontSize',12,'BoxStyle','full','Layer','top')
box on;hold on; axis tight;
% set(gca,'xTick',[ ],'yTick',[ ])
xlabel('$x$ [-]','FontSize',20,'interpreter','latex');
ylabel('$y$ [-]','FontSize',20,'interpreter','latex');
tit = strcat('Receptores:',' ',stringreceptores);
title(tit,'interpreter','latex','FontSize',20);
contourf(X,Y,rot90(MatDT,3),10), colormap('jet')
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = 'DT [-]';
c.Label.FontSize = 16;	
c.LineWidth = 1.2;
c.Limits = [round(min(min(MatDT)),1) round(max(max(MatDT)),1)];
filename2 = strcat(filename,'-',caso,caso2);
saveas(fig,fullfile(fpath,filename2),'epsc');
save(strcat('DT',caso2),'MatDT')