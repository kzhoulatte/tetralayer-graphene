clear all

kN =200;
thetaN =360;
kfilename = [num2str(kN),'x',num2str(thetaN),'_KP.mat']
load (kfilename);

tau  = 6.582e-13*2;
parfor_progress(100);

parfor i = 1:151
    parfor_progress;
    Eadd(1,i)=(i-1);
%     filename = ['inter_action/Eadd=',num2str(Eadd(1,i)),'_.mat'];
    filename = ['inter_action_200_360_Eadd=',num2str(Eadd(1,i)),'_.mat'];
    Rxx= [];
    
    Rxx = generate_Rxx_inter_1124(Kp,kr_all,filename,thetaN,tau);
    
    savefile = ['Rxx_inter_Rxx_1124_tau=',num2str(tau),'_inter_Eadd=',num2str(Eadd(1,i)),'_thetaN=',num2str(thetaN),'_.mat'];
    
    parsave(savefile, Rxx,'Rxx');
end

parfor_progress(0);