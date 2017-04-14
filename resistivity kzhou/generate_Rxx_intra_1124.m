function y=generate_Rxx_1124(Kp,kr_all,filename,thetaN,tau)


load(filename);

% load (kfilename);

cond_xx=zeros(size(Kp,1),301);

%  [y_val,y_vec,int_x]

gs =2;
gv =2;

% tau = 6.582e-13; %hbar/tau = 1meV

% KT = 25.7/298*1.5;

KT = 1;

dtheta = 2*pi/thetaN;

dk = kr_all(1,1);

for ief = 1:301
    Ef = (ief-151)/3;
    
    for kpp=1:size(Kp,1)
       for i = 1:size(int_x,1)
           if abs(y_val(i,i,kpp)-Ef)<= 5
               
               %% e^2/h * (gv*gs/2) *h*tau * (1/2pi)^2 * dthta * k* dk * vf^2 *(pfpe)
               
%                cond_xx(kpp,ief) = cond_xx(kpp,ief)+
%                kr_all(kpp,1)*int_x(i,i,kpp)^2*normpdf(y_val(i,i,kpp),Ef,2);
%                old simple form

%                

% %11_23                cond_xx(kpp,ief) = cond_xx(kpp,ief)+ (gv*gs/2)  *   tau    *   (2*pi)^-2 *dtheta*dk*kr_all(kpp,1)*abs(int_x(i,i,kpp))^2*(1/6.582e-16)*pfpe(y_val(i,i,kpp),Ef,KT);

                cond_xx(kpp,ief) = cond_xx(kpp,ief)+ (gv*gs/2)  *   tau    *   (2*pi)^-2 *dtheta*dk*kr_all(kpp,1)*abs(int_x(i,i,kpp))^2*(2*pi/(6.582e-16*1e3))*pfpe(y_val(i,i,kpp),Ef,KT); % consider the meV and eV
             
                
           end
       end
    end
end

y = cond_xx;

end
