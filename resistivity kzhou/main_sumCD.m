clear all

kN =200;
thetaN =360;
kfilename = [num2str(kN),'x',num2str(thetaN),'_KP.mat'];
load (kfilename);

for i = 1:151
    i
    Eadd(1,i)=(i-1);
%     filename = ['Rxx/Eadd=',num2str(Eadd(1,i)),'_.mat'];
    
%     filename = ['Rxx_1123/Rxx_1123_Eadd=',num2str(Eadd(1,i)),'_.mat'];
    filename = ['CD_CD_Eadd=',num2str(Eadd(1,i)),'_thetaN=',num2str(thetaN),'_.mat'];

    load (filename)
    
    CD_add(i,:) = sum(CD);
    
    clear CD;
end



surf(real(CD_add),'linestyle','none');

% Cxx_intra = Rxx_add;

save ('CD_1125.mat','CD_add');