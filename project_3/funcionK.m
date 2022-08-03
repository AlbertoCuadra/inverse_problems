function K = funcionK(x,y,L,T,N,D)
K = 0;% Inicializamos el valor de K para el sumatorio
for n=1:N
    K = K + (exp(-(n*pi)^2*D*T/L^2)*sin(n*pi*x/L)*sin(n*pi*y/L));
end
    K = 2/L*K;
end