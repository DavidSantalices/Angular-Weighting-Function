function s11_mie = getS11(m, r_part, lambda, ang)

% function S11=S11_Mie(m,r_part,lambda,ang)
%
% Función para calcular el parámetro S11 de scattering de una partícula
% para un ángulo de scattering dado.
%
% ARGUMENTOS DE ENTRADA:
% m: Índice de refracción complejo de la partícula
% r_part: Diámetro de la partícula (micras)
% lambda: Longitud de onda incidente en la partícula (micras) columna
% ang: Ángulo de scattering (rad) columna
%
% ARGUMENTOS DE SALIDA
% S11: Parámetro de scattering S11
%
% Autores: R. Magan y F. Betancourt     
% Date: 2017/07    Revision: 1.0 
% Copyright: LIR UC3M
% ________________________________________________________________________
% CÁLCULOS PREVIOS (Vector u, nmax y vector n)

x = 2*pi*r_part./lambda;
nmax=2.*ceil(x+4.*x.^(1/3)+2);
n=1:max(nmax);
nu=n+0.5;

% ________________________________________________________________________
% CÁLCULO DE LOS COEFICIENTES DE SCATTERING an y bn
n_arg = repmat(n,length(x),1); 
nu = repmat(nu,length(x),1);
mx = repmat(m.*x,1,max(nmax));
x_arg = repmat(x,1,max(nmax));

jmx=sqrt(pi./(2.*mx)).*besselj(nu,mx);
jx=sqrt(pi./(2.*x_arg)).*besselj(nu,x_arg);
yx=sqrt(pi./(2.*x_arg)).*bessely(nu,x_arg);
hx=jx+yx.*1i;

djx=[x.*sin(x)./x, x_arg(:,1:max(nmax)-1).*jx(:,1:max(nmax)-1)]-n_arg.*jx;
djmx=[m.*x.*sin(m.*x)./(m.*x), mx(:,1:max(nmax)-1).*jmx(:,...
    1:max(nmax)-1)]-n_arg.*jmx;
dhx=([x.*sin(x)./(x), x_arg(:,1:max(nmax)-1).*jx(:,1:max(nmax)-1)]+1i*...
    [-x.*cos(x)./x, x_arg(:,1:max(nmax)-1).*yx(:,1:nmax-1)])-n_arg.*hx;


%% hasta aqui guay las matrices
an=(m.^2.*jmx.*djx-jx.*djmx)./(m.^2.*jmx.*dhx-hx.*djmx);
bn=(jmx.*djx-jx.*djmx)./(jmx.*dhx-hx.*djmx);

% ________________________________________________________________________
% CÁLCULO DE LAS FUNCIONES DE AMPLITUD pn, tn, S1 y S2

s11_mie = zeros(length(x),length(ang));
M=(2.*n_arg+1)./(n_arg.*(n_arg+1));

for j=1:length(x)
    pn=ones(length(ang),nmax(j));
    tn=zeros(length(ang),nmax(j));
    S1=zeros(1,length(ang));%
    S2=zeros(1,length(ang));

    pn(:,2)=3.*cos(ang).*ones(length(ang),1);
    tn(:,1)=cos(ang).*ones(length(ang),1);
    tn(:,2)=3.*cos(2*ang).*ones(length(ang),1);
    
    for nn=3:nmax(j)
        pn(:,nn)=((2*nn-1).*cos(ang).*pn(:,nn-1)./(nn-1))-(nn.*...
            pn(:,nn-2)./(nn-1));
        tn(:,nn)=nn.*cos(ang).*pn(:,nn)-(nn+1).*pn(:,nn-1);
    end 
            
    for k =1:length(ang)
        S1(k)=sum(M(j,1:nmax(j)).*(an(j,1:nmax(j)).*pn(k,:)+...
            bn(j,1:nmax(j)).*tn(k,:)));
        S2(k)=sum(M(j,1:nmax(j)).*(an(j,1:nmax(j)).*tn(k,:)+...
            bn(j,1:nmax(j)).*pn(k,:)));
    end
    s11_mie(j,:)=(abs(S1).^2+abs(S2).^2)/2;
end

end