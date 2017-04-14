clear all 

kN =200;
thetaN =360;
kfilename = [num2str(kN),'x',num2str(thetaN),'_KP.mat']
load (kfilename);


parfor i = 1:151
    Eadd(1,i)=(i-1);
    
    i 
%     filename = ['inter_action/Eadd=',num2str(Eadd(1,i)),'_.mat'];
    filename = ['inter_action_200_360_Eadd=',num2str(Eadd(1,i)),'_.mat'];
    
    CD= [];
    CD = generate_charge_density(Kp,kr_all,filename,thetaN);
    
    savefile = ['CD_CD_Eadd=',num2str(Eadd(1,i)),'_thetaN=',num2str(thetaN),'_.mat'];
    
%     save (savefile,'Rxx');

    parsave_CD(savefile, CD , 'CD');

end