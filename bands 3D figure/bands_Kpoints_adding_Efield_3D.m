% Modified by Kuan Zhou April 12, 2017

clear all

Eadd=0;

delta1 = Eadd/2;     % H_tetra(1,1)
delta2 = Eadd/2;     % H_tetra(2,2)
delta3 = Eadd/6;     % H_tetra(3,3)
delta4 = Eadd/6;     % H_tetra(4,4)
delta5 = -1*Eadd/6;     % H_tetra(5,5)
delta6 = -1*Eadd/6;     % H_tetra(6,6)
delta7 = -1*Eadd/2;     % H_tetra(7,7)
delta8 = -1*Eadd/2;     % H_tetra(8,8)


% gamma0=3000;
% gamma1=390;
% gamma2=-16;
% gamma3=300;
% gamma4=40;
% gamma5=60;
% deltaAB=40.8;

% Supeng's parameters
gamma0=3000;
gamma1=400;
gamma2=-20/2;
gamma3=299;
gamma4=40;
gamma5=40/2;
deltaAB=50;

xi=1;

coolor='-k';

a=1.42*sqrt(3);
% 
G1=2*pi/a*[1,-1/sqrt(3)];
G2=2*pi/a*[0,2/sqrt(3)];

GA=(-G1*(2/3)-G2*(1/3))/40;
K2=(-G1*(1/3)+G2*(1/3))/40;
K1=[0,0,0];

 NK=500;
% 
 KX1=linspace(GA(1,1),K1(1,1),NK) ;
 KY1=linspace(GA(1,2),K1(1,2),NK) ;
% 
 KX2=linspace(K1(1,1),K2(1,1),NK) ;
 KY2=linspace(K1(1,2),K2(1,2),NK) ;

% Kp=[KX1,KX2;KY1,KY2];

kNx = NK;
kNy = NK;
% kfilename = [num2str(kNx),'x',num2str(kNy),'_KP.mat']
% load (kfilename);

Kx = linspace(-0.03,0.03,NK);
Ky = linspace(-0.03,0.03,NK);

% for kpp=1:size(Kp,2)
for i =1:kNx
    for j = 1:kNy

        %k=Kp(:,kpp)';
        k=[Kx(i),Ky(j)];

        V0PI=(sqrt(3)*a/2)*gamma0*(xi*k(1)+1i*k(2));

        V3PI=(gamma3/gamma0)*V0PI;
        V4PI=(gamma4/gamma0)*V0PI;


        % HG=[0,V0PI';V0PI,0];


        HG2up=[delta1,V0PI',-V4PI',V3PI;V0PI,deltaAB+delta2,gamma1,-V4PI';-V4PI,gamma1,deltaAB+delta3,V0PI';V3PI',-V4PI,V0PI,delta4];


        HG2down=[delta5,V0PI',-V4PI',V3PI;V0PI,deltaAB+delta6,gamma1,-V4PI';-V4PI,gamma1,deltaAB+delta7,V0PI';V3PI',-V4PI,V0PI,delta8];


        T2=[gamma2,0,0,0;0,gamma5,0,0;-V4PI,gamma1,gamma5,0;V3PI',-V4PI,0,gamma2];

        % T2=[gamma2,0,0,0;0,gamma5,0,0;-V4PI,gamma1,0,0;V3PI',-V4PI,0,0];

        % T2=[gamma2,0,0,0;0,gamma5,0,0;0,gamma1,gamma5,0;0,0,0,gamma2];
        % T2=zeros(4,4);

        HG4=[HG2up,T2;T2',HG2down];

        Egval(:,i,j)=sort(real(eig(HG4)));
    end
end


figure;

colormap([1 0 0;0.2 1 1]);

for i=1:4
    for p=1:kNx
        for q=1:kNy 
            val(p,q) = Egval(i,p,q);
        end
    end
    
    [kx,ky] = meshgrid(Kx,Ky);
    hSurf1 = surf(kx,ky,val,ones(NK)+1,'FaceAlpha',0.9);shading flat;
    hold on
end

for i=5:8
    for p=1:kNx
        for q=1:kNy 
            val(p,q) = Egval(i,p,q);
        end
    end
    [kx,ky] = meshgrid(Kx,Ky);
    hSurf2 = surf(kx,ky,val,ones(NK),'FaceAlpha',0.9);shading flat; 
    hold on
end


axis([-0.03 0.03 -0.03 0.03 -45 45]);

title('Band structure of tetralayer graphene');

view(-55, 17);

light('color',[0.2,0.7,1],'style','local','position',[1,0,0]);
light('color',[0.2,0.7,1],'style','local','position',[1,0,0]);
light('color',[0.2,0.7,1],'style','local','position',[0,1,0]);
light('color',[0.2,0.7,1],'style','local','position',[0,0,1]);
light('color',[0.2,0.7,1],'style','local','position',[0,1,0]);
lighting phong;

material shiny;
material([0.4 0 0.2 0.5 0.2]);

