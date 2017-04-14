clear all

kN =200;
thetaN =360;
kfilename = [num2str(kN),'x',num2str(thetaN),'_KP.mat']
% Load KP.mat file for kpoints

beep on;

parfor i = 1:151    % Eadd range
    i
    Eadd(1,i)=(i-1);
    y_val=[];
    y_vec=[];
    int_x=[];
    [y_val,y_vec,int_x] = generate_interaction(kfilename,Eadd(1,i));
    filename = ['inter_action_200_360_Eadd=',num2str(Eadd(1,i)),'_.mat'];
%     save(['inter_action_200_360_Ecenter/Eadd=',num2str(Eadd(1,i)),'_.mat'],'y_val','y_vec','int_x');

    parsave_interaction(filename, y_val,y_vec,int_x ,'y_val','y_vec','int_x');
    beep;

end

beep;