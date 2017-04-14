
Eadd = 1;
E_delta = 2;

delta1 = Eadd/2;    % H_tetra(1,1)
delta2 = Eadd/2;     % H_tetra(2,2)
delta3 = Eadd/6+E_delta;     % H_tetra(3,3)
delta4 = Eadd/6+E_delta;   % H_tetra(4,4)
delta5 = -1*Eadd/6+E_delta;    % H_tetra(5,5)
delta6 = -1*Eadd/6+E_delta;     % H_tetra(6,6)
delta7 = -1*Eadd/2;     % H_tetra(7,7)
delta8 = -1*Eadd/2;   % H_tetra(8,8)


gamma0=3000;
gamma1=390;
gamma2=-32/2;
gamma3=300;
gamma4=40;
gamma5=60;
deltaAB=40.8;

xi=1;

% coolor='-k'

a=1.42*sqrt(3);    % lattice constant of graphene

G1=2*pi/a*[1,-1/sqrt(3)];  % reciprocal lattice vectors (vector b1)
G2=2*pi/a*[0,2/sqrt(3)];  % vector b2


GA=(-G1*(2/3)-G2*(1/3))/40; % vector K  /40
K2=(-G1*(1/3)+G2*(1/3))/40; % vector Kp /40 
K1=[0,0,0]; 


load (kfilename);


for kpp=1:size(Kp,1)

k=Kp(kpp,:);
% fermi velocity times in-plane momentum; 
%Fermi velocity: (sqrt(3)*a/2)*gamma0
%In plane momentum: xi*k(1,1)+1i*k(1,2)

V0PI=(sqrt(3)*a/2)*gamma0*(xi*k(1,1)+1i*k(1,2));

V0PI_px = (sqrt(3)*a/2)*gamma0*(xi);

V3PI=(gamma3/gamma0)*V0PI;
V4PI=(gamma4/gamma0)*V0PI;

V3PI_px=(gamma3/gamma0)*V0PI_px;
V4PI_px=(gamma4/gamma0)*V0PI_px;

HG2up=[delta1,V0PI',-V4PI',V3PI;V0PI,deltaAB+delta2,gamma1,-V4PI';-V4PI,gamma1,deltaAB+delta3,V0PI';V3PI',-V4PI,V0PI,delta4];


HG2down=[delta5,V0PI',-V4PI',V3PI;V0PI,deltaAB+delta6,gamma1,-V4PI';-V4PI,gamma1,deltaAB+delta7,V0PI';V3PI',-V4PI,V0PI,delta8];


T2=[gamma2,0,0,0;0,gamma5,0,0;-V4PI,gamma1,gamma5,0;V3PI',-V4PI,0,gamma2];


HG4=[HG2up,T2;T2',HG2down]; % The total matrix 

HG2up_px=[0,V0PI_px,-V4PI_px,V3PI_px;V0PI_px,0,0,-V4PI_px;-V4PI_px,0,0,V0PI_px;V3PI_px,-V4PI_px,V0PI_px,0];

HG2down_px=[0,V0PI_px,-V4PI_px,V3PI_px;V0PI_px,0,0,-V4PI_px;-V4PI_px,0,0,V0PI_px;V3PI_px,-V4PI_px,V0PI_px,0];

T2_px = [0,0,0,0;0,0,0,0;-V4PI_px,0,0,0;V3PI_px,-V4PI_px,0,0];

HG4_px=[HG2up_px,T2_px;T2_px',HG2down_px];

V=[];
D=[];
interact_x=[];
interact_y=[];

[VV,DD] = eig(HG4);

%% only get the middle one

V = VV(:,3:6);

D = DD(3:6,3:6);


for ii=1:size(D,2)
    for jj=1:size(D,2)
     interact_x(ii,jj)=V(:,ii)'*HG4_px*V(:,jj);
    end
end

Egval(:,:,kpp)=real(D);
Egvec(:,:,kpp)=V;
int_x(:,:,kpp)=interact_x;

end

y=Egval;
y_vec=Egvec;
y_int_x=int_x;
