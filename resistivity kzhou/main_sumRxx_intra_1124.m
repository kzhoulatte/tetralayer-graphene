clear all

kN =200;
thetaN =360;
kfilename = [num2str(kN),'x',num2str(thetaN),'_KP.mat'];
load (kfilename);

tau  = 6.582e-13*2;

for i = 1:151
    i
    Eadd(1,i)=(i-1);
%     filename = ['Rxx/Eadd=',num2str(Eadd(1,i)),'_.mat'];
    
%     filename = ['Rxx_1123/Rxx_1123_Eadd=',num2str(Eadd(1,i)),'_.mat'];
    filename = ['Rxx_intra_Rxx_1124_tau=',num2str(tau),'_intra_Eadd=',num2str(Eadd(1,i)),'_thetaN=',num2str(thetaN),'_.mat'];

    load (filename)
    
    Rxx_add(i,:) = sum(Rxx);
    
    clear Rxx;
end

figure
surf(real(Rxx_add),'linestyle','none');

Cxx_intra = Rxx_add;

savefilename =['Rxx_intra_1125_360_tau=',num2str(tau),'.mat']

save (savefilename,'Cxx_intra');