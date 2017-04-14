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

CDall = [flipud(CD_add);CD_add];

% for i = 1:151
%     Eadd(1,i)=(i-1);
%     Eadd1(1,i)= (i-151);
% end

% Eadd_all=[Eadd1,Eadd];

% for i = 1:51
%     i;
%     Eadd(1,i)=(i-1);
% end

for i = 1:151
    Eadd(1,i)=(i-1);
    Eadd1(1,i)= (i-151);
end

Eadd =[Eadd1,Eadd];

for ief = 1:301
    Eff(1,ief) = (ief-151)/3;
end

colormap hot

save CDall;
save Eadd;
save Rall;  

 %figure
surf(CDall/1e12,Eadd,Rall,'Linestyle','none')

grid off
% colormap hot
% 
% 
% caxis([0,1])


% x = [0 0.001 0.005 0.01 0.05 1 2 3 4 5];
NNx =1000;

for i = 1:601
    nn(1,i)=(i+400)/1000;
end


for i = 1:size(nn,2)
    x(1,i) = (nn(1,i))^3;
end


Nx = length(x);

clim = [min(x) max(x)];

dx = min(diff(x));

y = clim(1):dx:clim(2);

for k=1:Nx-1
    y(y>x(k) & y<=x(k+1)) = x(k+1);
end

cmap = colormap(hot(Nx));

cmap2 = [...
interp1(x(:),cmap(:,1),y(:)) ...
interp1(x(:),cmap(:,2),y(:)) ...
interp1(x(:),cmap(:,3),y(:)) ...
];

colormap(cmap2)

caxis(clim)
c= colorbar;
c.Label.String = 'R_{xx} (k\Omega)';
az = 0;
el = 90;
view(az, el);


% figure 
% surf(CD_add,Eadd,abs(Cxx_intra),'Linestyle','none','EdgeAlpha','1')
% colormap hot
% 
% figure 
% 
% surf(CD_add,Eadd,abs(Cxx_inter),'Linestyle','none','EdgeAlpha','1')
% 

% figure
% plot(Call(1,:).^-1);
% hold on
% plot( Call(5,:).^-1);
% 
% plot(Call(10,:).^-1);
% 
% plot( Call(50,:).^-1);
% 
% hold off
% 
% 
% figure
% plot(Cxx_intra(1,:).^-1);
% hold on
% plot( Cxx_intra(5,:).^-1);
% 
% plot(Cxx_intra(10,:).^-1);
% 
% plot( Cxx_intra(50,:).^-1);
% 
% hold off


set(gcf, 'PaperPositionMode','auto')
    h_xlabel = get(gca, 'xlabel');
    h_ylabel = get(gca, 'ylabel');
%     set(gca,'linewidth',[4],'fontsize',[10],'position',[0.15,0.2,0.74,0.7])
      set(gca,'linewidth',[1],'fontsize',[12],'XTick',[-0.8   -0.4   0  0.4   0.8],'YTick',[-50 -25 0 25 50])
%       [-0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8]
    axis([-0.8 0.8 -150 150])
    x0=2;
    y0=1;
    width=5;
    height=4;
    set(gcf,'units','inches','position',[x0,y0,width,height])
    
    set(h_ylabel,'string','D (mV/nm)','fontsize',[16]);
    set(h_xlabel,'string','n(10^{12} cm^{-2})','fontsize',[16]);



beep

% save('data_test.mat','CDall','Eadd','Rall');