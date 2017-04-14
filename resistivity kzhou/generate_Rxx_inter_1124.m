function y=generate_Rxx_inter_1124(Kp,kr_all,filename,thetaN,tau)

load(filename);

% load (kfilename);

cond_xx_inter=zeros(size(Kp,1),301);

%  [y_val,y_vec,int_x]


gs =2;
gv =2;

% tau = 1e-12; %1ps

% tau = 6.582e-13; %hbar/tau = 1meV

% KT = 25.7/298;

KT = 1;

dtheta = 2*pi/thetaN;

dk = kr_all(1,1);


hbar_tau = 6.582e-16/tau*1e3; % meV


for ief = 1:301
    Ef = (ief-151)/3;
    
    for kpp=1:size(Kp,1)
       for i = 1:size(int_x,1)
           for j = 1:size(int_x,1)
           if (y_val(i,i,kpp)-Ef <= 5) && (y_val(i,i,kpp)-Ef > 0) && (Ef -y_val(j,j,kpp) <= 5) && (y_val(j,j,kpp)-Ef < 0) 
               
                    %% e^2/h * (gv*gs/2) *h* hbar * (1/2pi)^2 * dthta * k* dk * vf^2/(ei-ej) * (hbar/tau)/((ei-ej)^2+(hbar/tau)^2) 
               


%                

               cond_xx_inter(kpp,ief) = cond_xx_inter(kpp,ief)+  (gv*gs/2)  *  (2*pi) *   (2*pi)^-2    *dtheta*dk*kr_all(kpp,1)*int_x(i,j,kpp)*...
                   int_x(j,i,kpp)/(y_val(i,i,kpp)-y_val(j,j,kpp))*(hbar_tau)/((y_val(i,i,kpp)-y_val(j,j,kpp))^2+(hbar_tau)^2);
%                cond_xx(kpp,ief)+ (gv*gs/2)  *   tau    *   (2*pi)^-2 *dtheta*dk*kr_all(kpp,1)*abs(int_x(i,i,kpp))^2*(1/6.582e-16)*pfpe(y_val(i,i,kpp),Ef,KT);
               
               
               
%                cond_xx_inter(kpp,ief) = cond_xx_inter(kpp,ief)+ kr_all(kpp,1)*int_x(i,j,kpp)*int_x(j,i,kpp)*normpdf(y_val(i,i,kpp),y_val(j,j,kpp),2)/(y_val(i,i,kpp)-y_val(j,j,kpp));
           end
       end
    end
end

y = cond_xx_inter;
end
