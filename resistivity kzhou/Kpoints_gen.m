clear all

a=1.42*sqrt(3);

G1=2*pi/a*[1,-1/sqrt(3)];
G2=2*pi/a*[0,2/sqrt(3)];

% GA=(-G1*(2/3)-G2*(1/3))/20;
% K2=(-G1*(1/3)+G2*(1/3))/20;
% K1=[0,0,0];


% the reciprocal vecotr length 
GL = sqrt(sumsqr(G1));
KL = GL/40;

kN =200;
thetaN =360;

for i =1:(kN)
    kr  = i/kN*KL;
    for j = 1:thetaN
        angle = j/thetaN*360/180*pi;
        Kp(j+(i-1)*thetaN,:) = [kr *cos(angle),kr *sin(angle)];
        kr_all (j+(i-1)*thetaN,1) = kr;
    end
end

save([num2str(kN),'x',num2str(thetaN),'_KP.mat'],'Kp','kr_all');