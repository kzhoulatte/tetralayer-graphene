function y = generate_charge_density(Kp,kr_all,filename,thetaN)

load(filename);

CD =zeros(size(Kp,1),301);

dtheta = 2*pi/thetaN;

dk = kr_all(1,1);


for ief = 1:301
    Ef = (ief-151)/3;
        for kpp=1:size(Kp,1)
            for i = 1:size(int_x,1)
                if  y_val(i,i,kpp)<=Ef
                    CD(kpp,ief) = CD(kpp,ief)+ (2*pi)^-2 *dtheta*dk*kr_all(kpp,1);
                end
            end
        end


end

y =CD;

end