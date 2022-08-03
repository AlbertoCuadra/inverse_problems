function funcionK(x,y,T,N)
funcionK = 0;% Inicializamos el valor de K para el sumatorio
for n=1:N
    funcionK = funcionK + (exp(-(n*pi)^2*D*T/L^2)*sin(n*pi*x/L)*sin(n*pi*y/L));
end
    funcionK = 2/L*funcionK;
end