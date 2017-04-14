clear all

kN =200;
thetaN =360;
kfilename = [num2str(kN),'x',num2str(thetaN),'_KP.mat']
load (kfilename);

tau  = 6.582e-13*2;

for i = 1:151
    i
    Eadd(1,i)=(i-1);
    filename = ['Rxx_inter_Rxx_1124_tau=',num2str(tau),'_inter_Eadd=',num2str(Eadd(1,i)),'_thetaN=',num2str(thetaN),'_.mat'];

    load (filename)
    
    Rxx_add(i,:) = sum(Rxx);
    
    clear Rxx;
end

surf(real(Rxx_add),'linestyle','none');

Cxx_inter = Rxx_add;

savefilename =['Rxx_inter_1125_360_tau=',num2str(tau),'.mat']

save (savefilename,'Cxx_inter');

