clear all


tau  = 6.582e-13*2;

load(['Rxx_intra_1125_360_tau=',num2str(tau),'.mat']);

% load('Rxx_intra_1124_360.mat');
load(['Rxx_inter_1125_360_tau=',num2str(tau),'.mat']);

load ('CD_1125_mod.mat');

ratio = 2*pi;

Call =  abs(Cxx_intra +ratio*Cxx_inter)*3.869e-5*1e3;

R = Call.^-1;

Rall=[flipud(R);R];

CDall = [flipud(CD_add);CD_add]/1e12;

% for i = 1:151
%     Eadd(1,i)=(i-1);
%     Eadd1(1,i)= (i-151);
% end

% Eadd_all=[Eadd1,Eadd];

% for i = 1:51
%     i;
%     Eadd(1,i)=(i-1);
% end

for i = 1:51
    Eadd(1,i)=(i-1);
    Eadd1(1,i)= (i-51);
end

Eadd =[Eadd1,Eadd];

for ief = 1:301
    Eff(1,ief) = (ief-151)/3;
end

colormap hot

% figure
% surf(CDall/1e12,Eadd,Rall,'Linestyle','none')

grid off
% colormap hot
% 
% 
% caxis([0,1])



% surf(CD_add,Eadd,abs(Cxx_intra),'Linestyle','none','EdgeAlpha','1')
% colormap hot
% 
% figure 
% 
% surf(CD_add,Eadd,abs(Cxx_inter),'Linestyle','none','EdgeAlpha','1')
% 

%figure
plot(CDall(51,:),Rall(51,:));
hold on

plot(CDall(56,:),Rall(56,:));

plot(CDall(61,:),Rall(61,:));

plot(CDall(71,:),Rall(71,:));
% plot( Call(5,:).^-1);
% 
% plot(Call(10,:).^-1);
% 
% plot( Call(50,:).^-1);

hold off
% 
% 
figure
plot(CDall(51,:),Cxx_intra(1,:));
hold on


% hold on
% plot( Cxx_intra(5,:).^-1);
% 
% plot(Cxx_intra(10,:).^-1);
% 
% plot( Cxx_intra(50,:).^-1);
% 
hold off
figure
plot(CDall(51,:),Cxx_inter(1,:));

beep

% save('data_test.mat','CDall','Eadd','Rall');