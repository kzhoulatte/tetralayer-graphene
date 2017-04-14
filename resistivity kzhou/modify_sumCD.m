clear all

load ('CD_1125.mat');

for i = 1:151

    CD_add(i,:) = CD_add(i,:) - 0.0008699;
    
    clear CD;
    
end

CD_add =   CD_add*1e16;

save ('CD_1125_mod.mat','CD_add');

figure;

surf(real(CD_add),'linestyle','none');