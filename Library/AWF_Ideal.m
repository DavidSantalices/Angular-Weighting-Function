function [f,th] = AWF_Ideal(thL1,thL2,alL1,alL2,phL1,phL2,Nth,Nal,Nph,FoV,Source,D,Rang)
%Función que calcula la contribución de la geometría del problema cuando se da
%una función para el patrón de emisión de la fuente. El resultado dependerá
%por tanto del ángulo de scattering

syms alph theta phi

%Se expresan las coordenadas x, y, z en función de las nuevas coordenadas
%de nuestro problema (alpha, theta, phi)

x1=D/2 + (D*tand(alph)*(tand(alph)^2 - cosd(theta)^2 + 1)^(1/2) + D*tand(alph)^2*cosd(theta))/(2*tand(alph)*sind(theta)*(tand(alph)^2 + 1));
h1=(D*tand(alph)*(tand(alph)^2 - cosd(theta)^2 + 1)^(1/2) + D*tand(alph)^2*cosd(theta))/(2*sind(theta)*(tand(alph)^2 + 1));
y1=sind(phi).*h1;
z1=cosd(phi).*h1;

%Se hace un cambio de coordenadas para expresar las funciones de la fuente
%y los detectores en función de las nuevas coordenadas

FOV1=matlabFunction(FoV(x1,y1,z1));
% [xx,yy]=meshgrid([-100:100],[-100:100]);
% 
% 
% for i=9.6:10:100
% A=FoV(xx,yy,i).*Source(xx,yy,i);
% figure
% pcolor(A)
% end


%El objetivo de esta aproximación es realizar un tratamiento continuo del
%problema y por lo tanto habrá que realizar una integral. Al cambiar las
%coordenadas por otras cuyos diferenciales no son uniformes a lo largo del
%espacio es necesario calcular el jacobiano
v1=[x1,y1,z1];
u1=[theta alph phi];
H1 = simplify(det(jacobian (v1,u1)));
J1=matlabFunction(H1);

xx1= matlabFunction(x1);
zz1= matlabFunction(z1);
yy1= matlabFunction(y1);

%Se ha divido el espacio en dos zonas una a la derecha del punto medio que
%une fuente y detector y otra a la izquierda de ese punto medio. En este
%caso se realiza los mismos pasos que anteriormente pero para la otra zona
%del espacio
x2=D/2 - (D*tand(alph)*(tand(alph)^2 - cosd(theta)^2 + 1)^(1/2) - D*tand(alph)^2*cosd(theta))/(2*tand(alph)*sind(theta)*(tand(alph)^2 + 1));
h2=-(D*tand(alph)*(tand(alph)^2 - cosd(theta)^2 + 1)^(1/2) - D*tand(alph)^2*cosd(theta))/(2*sind(theta)*(tand(alph)^2 + 1));
y2=sind(phi).*h2;
z2=cosd(phi).*h2;


FOV2=matlabFunction(FoV(x2,y2,z2));

v2=[x2,y2,z2];
u2=[theta alph phi];
H2 = simplify(det(jacobian (v2,u2)));
J2=matlabFunction(H2);

xx2= matlabFunction(x2);
zz2= matlabFunction(z2);
yy2= matlabFunction(y2);


%Se crea un meshgrid para calcular la contribución de la geometría al
%problema en cada punto del espacio. Como se ha mencionado anteriormente se
%ha dividido el espacio en dos zonas, una dada por angulos alpha que van de
%0 a 90 y otra dada por angulos alpha de 0 a -90
[ALP1,PH1,TH1]=meshgrid(linspace(0,alL2,Nal),linspace(phL1,phL2,Nph),linspace((180-thL1),(180-thL2),Nth));
[ALP2,PH2,TH2]=meshgrid(linspace(0,-alL2,Nal),linspace(phL1,phL2,Nph),linspace((180-thL1),(180-thL2),Nth));

X1=xx1(ALP1,TH1);
Y1=yy1(ALP1,PH1,TH1);
Z1=zz1(ALP1,PH1,TH1);

X2=xx2(ALP2,TH2);
Y2=yy2(ALP2,PH2,TH2);
Z2=zz2(ALP2,PH2,TH2);

RA1=matlabFunction(Rang(x1,y1,z1));
RA2=matlabFunction(Rang(x2,y2,z2));

Fsource1=matlabFunction(Source(x1,y1,z1));
Fsource2=matlabFunction(Source(x2,y2,z2));

r2_1=1./((xx1(ALP1,TH1)-D).^2+yy1(ALP1,PH1,TH1).^2+zz1(ALP1,PH1,TH1).^2);
r2_2=1./((xx2(ALP2,TH2)-D).^2+yy2(ALP2,PH2,TH2).^2+zz2(ALP2,PH2,TH2).^2);

F1=J1(ALP1,TH1).*Fsource1(ALP1,PH1,TH1).*FOV1(ALP1,PH1,TH1).*RA1(ALP1,PH1,TH1).*r2_1;
F2=J2(ALP2,TH2).*Fsource2(ALP2,PH2,TH2).*FOV2(ALP2,PH2,TH2).*RA2(ALP1,PH1,TH1).*r2_2;


ALP=cat(2,ALP1,ALP2);
PH=cat(2,PH1,PH2);
TH=cat(2,TH1,TH2);

% X=cat(2,X1,X2);
% Y=cat(2,Y1,Y2);
% Z=cat(2,Z1,Z2);

% [X,Y,Z]=meshgrid(linspace(0,D),linspace(-D,D),linspace(0,D));
% 
% FF=Source(X,Y,Z).*Rang(X,Y,Z).*FoV(X,Y,Z).*(1./((X-D).^2+Y.^2+Z.^2));
% % 
% 
% [al1,ph1]=meshgrid(linspace(0,90),linspace(0,90));
% [al2,ph2]=meshgrid(linspace(-90,-1),linspace(0,90));
% 
% figure
% hold on
% for t=85:10:160
% xslice=cat(2,xx1(al1,t),xx2(al2,t));
% yslice=cat(2,yy1(al1,ph1,t),yy2(al2,ph2,t));
% zslice=cat(2,zz1(al1,ph1,t),zz2(al2,ph2,t));
% % [xslice,yslice] = meshgrid(linspace(0,D,200),linspace(0,D,200));
% % zslice = 0*yslice;
% 
% h=slice(X,Y,Z,FF,xslice,yslice,zslice);
% h.EdgeColor='None'
% % axis equal
% % hold on
% % h=slice(X,Y,Z,FF,[],0,[])
% h.FaceAlpha=0.8;
% end
% %Calculamos la integral en ambas zonas y sumamos ambas contribuciones
% axis equal
% 
% xlabel('x (mm)')
% ylabel('y (mm)')
% zlabel('z (mm)')


F1(isnan(F1))=0;

F2(isnan(F2))=0;

f1=squeeze(trapz(trapz(F1,2),1))*abs((ALP1(2,2,2)-ALP1(1,1,1)).*(PH1(2,2,2)-PH1(1,1,1)));

f2=squeeze(trapz(trapz(F2,2),1))*abs((ALP2(2,2,2)-ALP2(1,1,1)).*(PH2(2,2,2)-PH2(1,1,1)));

f=f1+f2;

% figure
% pcolor(squeeze(X1()),squeeze(Y1()),squeeze())
% hold on
% pcolor()

[a,b,c]=size(TH1);

%Cálculo del volumen de interacción
%Como la interacción es "continua" en todo el volumen se debe establecer un
%umbral, en este caso se va a tomar el volumen donde la interacción tiene una
%intensidad mayor del 20% del máximo.
%


th=linspace(thL1,thL2,Nth);
end

