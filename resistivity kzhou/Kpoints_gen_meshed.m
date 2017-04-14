clear all

a=1.42*sqrt(3);

G1=2*pi/a*[1,-1/sqrt(3)];
G2=2*pi/a*[0,2/sqrt(3)];

% GA=(-G1*(2/3)-G2*(1/3))/20;
% K2=(-G1*(1/3)+G2*(1/3))/20;
% K1=[0,0,0];

% the reciprocal vecotr length 

NK = 100;
kNx = 100;
kNy = 100; 

Kx = linspace(-1,1,NK);
Ky = linspace(-1,1,NK);

for i =1:(kNx)
    for j = 1:kNy
        Kp(j+(i-1)*100,:) = [Kx(i),Ky(j)];
    end
end

save([num2str(kNx),'x',num2str(kNy),'_KP.mat'],'Kp');