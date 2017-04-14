
%**************************
%   Code to calculate spectrum of 4LG_Kp
%   with gamma1 gamma3, gamma4
%   gamma2 and gamma5 involved.
%
%***************************

%*********************************
% choosen base eign wave
%   |phi_1>=(|0>,0,0,0,0,0,0,0);
%   |phi_2>=(0,0,|0>,0,0,0,0,0);
%   |phi_3n>=1/sqrt(2)(|n+1>,|n>,0,0,0,0,0,0);
%   |phi_4n>=1/sqrt(2)(0,0,|n+1>,|n>,0,0,0,0);
%   |phi_5n>=1/sqrt(2)(-|n+1>,|n>,0,0,0,0,0,0);
%   |phi_6n>=1/sqrt(2)(0,0,-|n+1>,|n>,0,0,0,0);
%
%   |phi_7>=(0,0,0,0,|0>,0,0,0);
%   |phi_8>=(0,0,0,0,0,0,|0>,0);
%   |phi_9n>=1/sqrt(2)(0,0,0,0,|n+1>,|n>,0,0);
%   |phi_10n>=1/sqrt(2)(0,0,0,0,0,0,|n+1>|n>,);
%   |phi_11n>=1/sqrt(2)(0,0,0,0,-|n+1>,|n>,0,0);
%   |phi_12n>=1/sqrt(2)(0,0,0,0,0,0,-|n+1>,|n>);
%**********************************

clear all;

Emin = -80;
Emax = 80;
N = 50; % cutoff Landau Level
Bmin = 0.01;
Bmax = 12; % Maximum B field
stepsB=501;     % steps of B field for calculation

%----------------------------------
%
%    hopping parameters for calculations
%----------------------------------
gamma0 = 3000; 
r1BLG = 390;
r1Inter = 320;
%gamma1 = 390; % gamma1 value in the unit of meV
%gamma1 = 0; % gamma1 value in the unit of meV

gamma3 =300;
gamma4 = 40;
% gamma3 =0;
% gamma4 = 0;

gamma2 = -10;
gamma5 = 20;
% gamma2 = 0;
% gamma5 = 0;

deltaAB = 40.8; % energy difference between stacked and non-stacked atoms
delta_middle =0;    % potential difference between the middle two layers and the outer layers
Eadd=45; % The electric field applied to TLG

%**************************************
%     Electric field added
%
%**************************************

delta1 = Eadd/2;    % H_tetra(1,1)
delta2 = Eadd/2;     % H_tetra(2,2)
delta3 = Eadd/6;     % H_tetra(3,3)
delta4 = Eadd/6;   % H_tetra(4,4)
delta5 = -1*Eadd/6;    % H_tetra(5,5)
delta6 = -1*Eadd/6;     % H_tetra(6,6)
delta7 = -1*Eadd/2;     % H_tetra(7,7)
delta8 = -1*Eadd/2;   % H_tetra(8,8)

%**************************************
%    The middle layer energy difference
%    whe E=0
%
%***************************************

delta3 =  delta3 + delta_middle;
delta4 =  delta4 + delta_middle;
delta5 =  delta5 + delta_middle;
delta6 =  delta6 + delta_middle;


%**************************************
%    Energy difference between AB sublattices
%
%
%**************************************

delta2 = delta2+deltaAB;
delta3 = delta3+deltaAB;
delta6 = delta6+deltaAB;
delta7 = delta7+deltaAB;



%
%        H,    H_i
% Ht = (           )   
%        H_i,   H

E = zeros(stepsB,13+8*N);

%-----------------
%   Index for convenience
%----------------
index1 = 1;
index2 = 2;

min3 = 3;
max3 = 3+N;

min4 = 4+N;
max4 = 4+2*N;

min5 = 5+2*N;
max5 = 5+3*N;

min6 = 6+3*N;
max6 = 6+4*N;

index7 = 7+4*N;
index8 = 8+4*N;

min9 = 9+4*N;
max9 = 9+5*N;

min10 = 10+5*N;
max10 = 10+6*N;

min11 = 11+6*N;
max11 = 11+7*N;

min12 = 12+7*N;
max12 = 12+8*N;


 for step = 1:stepsB
    
    gamma1 = r1BLG;
    
    B = (step-1)/(stepsB-1)*(Bmax-Bmin) + Bmin;

    H = zeros(6+4*N);  % Hamiltonian of BLG part
    Ht = zeros(12+8*N);
    
    
    
%     epslon0 = 35.223*sqrt(abs(B));  % Epsilon0=36.2*sqrt(B) in the unit of meV
    epslon0 = 36.26*sqrt(abs(B));
%*******************************************************
% Bilayer blocks
    
    % <phi_1|H|phi>
    m=index1;
    LandauM = m-index1;

    %   <phi_1|H|phi_1> == 0
    for n = index1:index1
        LandauN = n-index1;
        H(m,n) = H(m,n) + delta1*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end                
    %   <phi_1|H|phi_2> == 0

    %  <phi_1|H|phi_3_n> == 0

    %  <phi_1|H|phi_4_n>
    for n = min4:max4
        LandauN = n-min4;
        H(m,n) = H(m,n) - 1/sqrt(2)*(gamma3/gamma0)*epslon0*sqrt(LandauN)*delta(LandauN-1, LandauM);
        H(n,m) = H(m,n);
    end


    % <phi_1|H|phi_5_n> ==0

    % <phi_1|H|phi_6_n>
    for n = min6:max6
        LandauN = n-min6;
        H(m,n) = H(m,n) - 1/sqrt(2)*(gamma3/gamma0)*epslon0*sqrt(LandauN)*delta(LandauN-1, LandauM);
        H(n,m) = H(m,n);
    end


% <phi_2|H|phi>
    m=index2;
    LandauM = m-index2;

    % <phi_2|H|phi_2> == 0
    for n = index2:index2
        LandauN = n-index2;
        H(m,n) = H(m,n) + delta3*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end    

    %  <phi_2|H|phi_3_n>
    for n = min3:max3
        LandauN = n-min3;
        H(m,n) = H(m,n) + 1/sqrt(2)*(gamma4/gamma0*epslon0*sqrt(LandauN+1)+gamma1)*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end

    %<phi_2|H|phi_4_n> ==0


    % <phi_2|H|phi_5_n>
    for n = min5:max5
        LandauN = n-min5;
        H(m,n) = H(m,n) + 1/sqrt(2)*((-1)*gamma4/gamma0*epslon0*sqrt(LandauN+1)+gamma1)*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end

    % <phi_2|H|phi_6_n> ==0

                
     
    %<phi_3_m|H|phi>
    for m = min3:max3
        
        LandauM = m-min3;
        % <phi_3_m|H|phi_3_n>
        for n = m:max3
            LandauN = n-min3;
            H(m,n) = H(m,n) + ((delta1+delta2)/2-epslon0*sqrt(LandauN+1))*delta(LandauM,LandauN);
            H(n,m) = H(m,n);
        
        end


        % <phi_3_m|H|phi_4_n>
       for n = min4:max4
            LandauN = n-min4;
            H(m,n) = H(m,n) + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+1)...
                        -gamma3/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN-1)...
                         +(gamma1+gamma4/gamma0*epslon0*sqrt(LandauN+1))*delta(LandauM,LandauN+1));
            H(n,m) = H(m,n);
        end


        % <phi_3_m|H|phi_5_n> == 0
        for n = min5:max5
            LandauN = n-min5;
            H(m,n) = H(m,n) + 1/2*(delta2*delta(LandauM, LandauN)-delta1*delta(LandauM+1, LandauN+1));
            H(n,m) = H(m,n);
        end

        % <phi_3_m|H|phi_6_n>
        for n = min6:max6
            LandauN = n-min6;
            H(m,n) = H(m,n) + 1/2*(-gamma4/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+1)...
                            -gamma3/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN-1)...
                            +(-gamma1+gamma4/gamma0*epslon0*sqrt(LandauN+1))*delta(LandauM, LandauN+1));
            H(n,m) = H(m,n);

        end


    end;


    %<phi_4_m|H|phi>
    for m = min4:max4
        
        LandauM = m-min4;
        % <phi_4_m|H|phi_4_n>
        for n = m:max4
            LandauN = n-min4;
            H(m,n) = H(m,n) - epslon0*sqrt(LandauN+1)*delta(LandauM, LandauN)+1/2*(delta3+delta4)*delta(LandauM, LandauN);
            H(n,m) = H(m,n);
        end

        % <phi_4_m|H|phi_5_n>
        for n = min5:max5
            LandauN = n-min5;
            H(m,n) = H(m,n) + 1/2*((-gamma4/gamma0*epslon0*sqrt(LandauN+1)+gamma1)*delta(LandauM+1, LandauN)...
                               + gamma3/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+2)...
                               + gamma4/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN));
            H(n,m) = H(m,n);

        end
        % <phi_4_m|H|phi_6_n> == 0 
            for n = min6:max6
            LandauN = n-min6;
            H(m,n) = H(m,n) + 1/2*(delta4*delta(LandauM, LandauN)...
                               -delta3*delta(LandauM+1, LandauN+1));
            H(n,m) = H(m,n);
            end
        

    end;

    % <phi_5_m|H|phi>
    for m = min5:max5
        
        LandauM = m-min5;
        % <phi_5_m|H|phi_5_n>
        for n = m:max5
            LandauN = n-min5;
            H(m,n) = H(m,n) + 1/2*(+2*epslon0*sqrt(LandauN+1)+delta1+delta2)*delta(LandauM, LandauN);
            H(n,m) = H(m,n);

        end

        % <phi_5_m|H|phi_6_n>
        for n = min6:max6
            LandauN = n-min6;
            H(m,n) = H(m,n)  + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM, LandauN+1)...
                            +gamma3/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN-1)...
                            +(-gamma1+gamma4/gamma0*epslon0*sqrt(LandauN+1))*delta(LandauM, LandauN+1));
            H(n,m) = H(m,n);

        end        
    end;

    % <phi_6_m|H|phi>
    for m = min6:max6
        
        LandauM = m-min6;
        % <phi_6_m|H|phi_6_n>
        for n = m:max6
            LandauN = n-min6;
            H(m,n) =H(m,n) + (epslon0*sqrt(LandauN+1)+1/2*(delta3+delta4))*delta(LandauM, LandauN);
            H(n,m) = H(m,n);

        end

    end
    
    % Filling the two BLG blocks
    Ht(1:6+4*N,1:6+4*N)=H;
    
H = zeros(6+4*N);    



%*******************************************************
%another Bilayer blocks


    % <phi_1|H|phi>
    m=index1;
    LandauM = m-index1;

    %   <phi_1|H|phi_1> == 0
    for n = index1:index1
        LandauN = n-index1;
        H(m,n) = H(m,n) + delta5*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end                
    %   <phi_1|H|phi_2> == 0

    %  <phi_1|H|phi_3_n> == 0

    %  <phi_1|H|phi_4_n>
    for n = min4:max4
        LandauN = n-min4;
        H(m,n) = H(m,n) - 1/sqrt(2)*(gamma3/gamma0)*epslon0*sqrt(LandauN)*delta(LandauN-1, LandauM);
        H(n,m) = H(m,n);
    end


    % <phi_1|H|phi_5_n> ==0

    % <phi_1|H|phi_6_n>
    for n = min6:max6
        LandauN = n-min6;
        H(m,n) = H(m,n) - 1/sqrt(2)*(gamma3/gamma0)*epslon0*sqrt(LandauN)*delta(LandauN-1, LandauM);
        H(n,m) = H(m,n);
    end


% <phi_2|H|phi>
    m=index2;
    LandauM = m-index2;

    % <phi_2|H|phi_2> == 0
    for n = index2:index2
        LandauN = n-index2;
        H(m,n) = H(m,n) + delta7*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end    

    %  <phi_2|H|phi_3_n>
    for n = min3:max3
        LandauN = n-min3;
        H(m,n) = H(m,n) + 1/sqrt(2)*(gamma4/gamma0*epslon0*sqrt(LandauN+1)+gamma1)*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end

    %<phi_2|H|phi_4_n> ==0


    % <phi_2|H|phi_5_n>
    for n = min5:max5
        LandauN = n-min5;
        H(m,n) = H(m,n) + 1/sqrt(2)*((-1)*gamma4/gamma0*epslon0*sqrt(LandauN+1)+gamma1)*delta(LandauN, LandauM);
        H(n,m) = H(m,n);
    end

    % <phi_2|H|phi_6_n> ==0

                
     
    %<phi_3_m|H|phi>
    for m = min3:max3
        
        LandauM = m-min3;
        % <phi_3_m|H|phi_3_n>
        for n = m:max3
            LandauN = n-min3;
            H(m,n) = H(m,n) + ((delta5+delta6)/2-epslon0*sqrt(LandauN+1))*delta(LandauM,LandauN);
            H(n,m) = H(m,n);
        
        end


        % <phi_3_m|H|phi_4_n>
       for n = min4:max4
            LandauN = n-min4;
            H(m,n) = H(m,n) + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+1)...
                        -gamma3/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN-1)...
                         +(gamma1+gamma4/gamma0*epslon0*sqrt(LandauN+1))*delta(LandauM, LandauN+1));
            H(n,m) = H(m,n);
        end


        % <phi_3_m|H|phi_5_n> == 0
        for n = min5:max5
            LandauN = n-min5;
            H(m,n) = H(m,n) + 1/2*(delta6*delta(LandauM, LandauN)-delta5*delta(LandauM+1, LandauN+1));
            H(n,m) = H(m,n);
        end

        % <phi_3_m|H|phi_6_n>
        for n = min6:max6
            LandauN = n-min6;
            H(m,n) = H(m,n) + 1/2*(-gamma4/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+1)...
                            -gamma3/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN-1)...
                            +(-gamma1+gamma4/gamma0*epslon0*sqrt(LandauN+1))*delta(LandauM, LandauN+1));
            H(n,m) = H(m,n);

        end


    end;


    %<phi_4_m|H|phi>
    for m = min4:max4
        
        LandauM = m-min4;
        % <phi_4_m|H|phi_4_n>
        for n = m:max4
            LandauN = n-min4;
            H(m,n) = H(m,n) - epslon0*sqrt(LandauN+1)*delta(LandauM, LandauN)+1/2*(delta7+delta8)*delta(LandauM, LandauN);
            H(n,m) = H(m,n);
        end

        % <phi_4_m|H|phi_5_n>
        for n = min5:max5
            LandauN = n-min5;
            H(m,n) = H(m,n) + 1/2*((-gamma4/gamma0*epslon0*sqrt(LandauN+1)+gamma1)*delta(LandauM+1, LandauN)...
                               + gamma3/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+2)...
                               + gamma4/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN));
            H(n,m) = H(m,n);

        end
        % <phi_4_m|H|phi_6_n> == 0 
            for n = min6:max6
            LandauN = n-min6;
            H(m,n) = H(m,n) + 1/2*(delta8*delta(LandauM, LandauN)...
                               -delta7*delta(LandauM+1, LandauN+1));
            H(n,m) = H(m,n);
            end
        

    end;

    % <phi_5_m|H|phi>
    for m = min5:max5
        
        LandauM = m-min5;
        % <phi_5_m|H|phi_5_n>
        for n = m:max5
            LandauN = n-min5;
            H(m,n) = H(m,n) + 1/2*(+2*epslon0*sqrt(LandauN+1)+delta5+delta6)*delta(LandauM, LandauN);
            H(n,m) = H(m,n);

        end

        % <phi_5_m|H|phi_6_n>
        for n = min6:max6
            LandauN = n-min6;
            H(m,n) = H(m,n)  + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM, LandauN+1)...
                            +gamma3/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1, LandauN-1)...
                            +(-gamma1+gamma4/gamma0*epslon0*sqrt(LandauN+1))*delta(LandauM, LandauN+1));
            H(n,m) = H(m,n);

        end        
    end;

    % <phi_6_m|H|phi>
    for m = min6:max6
        
        LandauM = m-min6;
        % <phi_6_m|H|phi_6_n>
        for n = m:max6
            LandauN = n-min6;
            H(m,n) =H(m,n) + (epslon0*sqrt(LandauN+1)+1/2*(delta7+delta8))*delta(LandauM, LandauN);
            H(n,m) = H(m,n);

        end

    end




%%%%%%
Ht(7+4*N:12+8*N,7+4*N:12+8*N)=H;    
    
    
%*******************************************************   
%   Interaction part of the two BLG

%  with gamma1 between A2 site and B3 site
%   
%       
    gamma1 = r1Inter;
    % <phi_1|H_i|phi> == 0
    
    % <phi_2|H_i|phi> ~= 0
    
    %<phi_2|H_i|phi>
        for m=index2;
        LandauM = m-index2;
        
    %<phi_2|H_i|phi9n>    
         for n = min9:max9
             LandauN = n-min9;        
             Ht(m,n) = Ht(m,n) + 1/sqrt(2)*gamma1*delta(LandauM, LandauN);
             Ht(n,m) = Ht(m,n);
         end            
 
    %<phi_2|H_i|phi11n>    
         for n = min11:max11
             LandauN = n-min11;        
             Ht(m,n) = Ht(m,n) + 1/sqrt(2)*gamma1*delta(LandauM, LandauN);
             Ht(n,m) = Ht(m,n);
         end            
        end
        
    % <phi_3n|H_i|phi> == 0
    
    % <phi_4n|H_i|phi>
        for m = min4: max4
        LandauM = m-min4;
        % <phi_4m|H_i|phi_4n> == 0
        
        % <phi_4m|H_i|phi_5n> == 0
        
        % <phi_4m|H_i|phi_6n> == 0
        
        % <phi_4m|H_i|phi_7n>==0
        
         % <phi_4m|H_i|phi_8n> == 0;
         
         % <phi_4m|H_i|phi_9n> == 0;
         for n = min9:max9
             LandauN = n-min9;
             Ht(m,n) = Ht(m,n) + 1/2*gamma1*delta(LandauM+1, LandauN);
             Ht(n,m) = Ht(m,n);
         end    
         
         % <phi_4m|H_i|phi_10n> == 0;
         
         % <phi_4m|H_i|phi_11n> == 0;
         for n = min11:max11
             LandauN = n-min11;
             Ht(m,n) =Ht(m,n) +  1/2*gamma1*delta(LandauM+1, LandauN);
             Ht(n,m) = Ht(m,n);
         end
         % <phi_4m|H_i|phi_12n> == 0;         
        end
        
        
    % <phi_5n|H_i|phi> == 0
    
    % <phi_6n|H_i|phi> 
     for m = min6: max6
        LandauM = m-min6;   
        % <phi_6m|H_i|phi_7n>==0

        % <phi_6m|H_i|phi_8n> == 0;
         
         % <phi_6m|H_i|phi_9n> == 0;
         for n = min9:max9
             LandauN = n-min9;
             Ht(m,n) = Ht(m,n) +  -1/2*gamma1*delta(LandauM+1, LandauN);
             Ht(n,m) = Ht(m,n);
         end
         % <phi_6m|H_i|phi_10n> == 0;
         
         % <phi_6m|H_i|phi_11n> == 0;
         for n =  min11:max11
             LandauN = n-min11;
             Ht(m,n) = Ht(m,n) +  -1/2*gamma1*delta(LandauM+1, LandauN);
             Ht(n,m) = Ht(m,n);
         end
         
         % <phi_6m|H_i|phi_12n> == 0;
     end
    % <phi_7n|H_i|phi> == 0
    % <phi_8n|H_i|phi> == 0
    % <phi_9n|H_i|phi> == 0
    % <phi_10n|H_i|phi> == 0
    % <phi_11n|H_i|phi> == 0
    % <phi_12n|H_i|phi> == 0
    
    
%***********************************************************
%***********************************************************
%
%       coupling between A2-A3 and B2-B3
%
%
%
%************************************************************
        % <phi_1|H|phi> ==0
        % <phi_2|H|phi>
        m=index2;
        LandauM = m-index2;
        % <phi_2|H|phi_9_n>
        for n = min9:max9
            LandauN = n-min9;
            Ht(m,n) = Ht(m,n) + 1/sqrt(2)*gamma4/gamma0*epslon0*sqrt(LandauN+1)*delta(LandauN,LandauM);
            Ht(n,m) = Ht(m,n);
        end
        % <phi_2|H|phi_11_n>
        for n = min11:max11
            LandauN = n-min11;
            Ht(m,n) = Ht(m,n) - 1/sqrt(2)*gamma4/gamma0*epslon0*sqrt(LandauN+1)*delta(LandauN,LandauM);
            Ht(n,m) = Ht(m,n);
        end
        
    % <phi_4_m|H|phi>
    for m = min4:max4
        LandauM = m-min4;
   
         % <phi_4_m|H|phi_7_n>
       n = index7;
            LandauN = n-index7;
            Ht(m,n) = Ht(m,n) - 1/sqrt(2)*gamma3/gamma0*epslon0*delta(LandauN+1,LandauM);
            Ht(n,m) = Ht(m,n);
       
         
         
        
        % <phi_4_m|H|phi_9_n>
        for n = min9:max9
            LandauN = n-min9;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN+1)*delta(LandauM+1,LandauN)...
                                        -gamma3/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+2)...
                                        +gamma4/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1,LandauN));
            Ht(n,m) = Ht(m,n);
        end
        % <phi_4|H|phi_11_n>
        for n = min11:max11
            LandauN = n-min11;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1,LandauN)...
                                        +gamma3/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+2)...
                                        -gamma4/gamma0*epslon0*sqrt(LandauN+1)*delta(LandauM+1,LandauN));
            Ht(n,m) = Ht(m,n);
        end
    end
        
        
     % <phi_6_m|H|phi>
    for m = min6:max6
        LandauM = m-min6;
        
       % <phi_6_m|H|phi_7>
        n=index7;
        LandauN = n-index7;
        Ht(m,n) = Ht(m,n) - 1/sqrt(2)*gamma3/gamma0*epslon0*delta(LandauM,LandauN+1);
        Ht(n,m) = Ht(m,n);  
        
        
        
        % <phi_6_m|H|phi_9_n>
        for n = min9:max9
            LandauN = n-min9;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN)*delta(LandauM+1,LandauN)...
                                        -gamma3/gamma0*epslon0*sqrt(LandauN+2)*delta(LandauM,LandauN+2)...
                                        -gamma4/gamma0*epslon0*sqrt(LandauN+1)*delta(LandauM+1,LandauN));
            Ht(n,m) = Ht(m,n);
        end
        % <phi_6|H|phi_11_n>
        for n = min11:max11
            LandauN = n-min11;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma4/gamma0*epslon0*sqrt(LandauN)*delta(LandauM,LandauN-1)...
                                        +gamma3/gamma0*epslon0*sqrt(LandauN+1)*delta(LandauM,LandauN+2)...
                                        +gamma4/gamma0*epslon0*sqrt(LandauN+1)*delta(LandauM+1,LandauN));
            Ht(n,m) = Ht(m,n);
        end
         
    end
    
    
%***************************************************************************
%
%
%
%   hopping parameter gamma2 and gamma5 involved.
%
%
%
%
%
%*****************************************************************************

    %<phi_1|H|phi> is zero except <phi_1|H|phi_7>
        m = index1;
        LandauM = m-index1;
        n=  index7;
        LandauN = n-index7;
        Ht(m,n) = Ht(m,n) + gamma2;
        Ht(n,m) = Ht(m,n);
        
    %<phi_2|H|phi> is zero except <phi_1|H|phi_8>
        m = index2;
        LandauM = m-index2;
        n=  index8;
        LandauN = n-index8;
        Ht(m,n) = Ht(m,n) + gamma5;
        Ht(n,m) = Ht(m,n);

    %<phi_3_m|H|phi> is zero except <phi_3_m|H|phi_9n> and <phi_3_m|H|phi_11n>
    
    for m = min3:max3
        LandauM = m-min3;
        % <phi_3_m|H|phi_9n>
        for n = min9:max9
            LandauN = n-min9;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma2+gamma5)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
        % <phi_3_m|H|phi_11n>
        for n = min11:max11
            LandauN = n-min11;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma5-gamma2)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
    end
    
    %<phi_4_m|H|phi> is zero except <phi_4_m|H|phi_10n> and <phi_4_m|H|phi_12n>
    
    for m = min4:max4
        LandauM = m-min4;
        % <phi_4_m|H|phi_9n>
        for n = min10:max10
            LandauN = n-min10;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma2+gamma5)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
        % <phi_4_m|H|phi_12n>
        for n = min12:max12
            LandauN = n-min12;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma2-gamma5)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
    end
        
    %<phi_5_m|H|phi> is zero except <phi_5_m|H|phi_9n> and <phi_5_m|H|phi_11n>
    
    for m = min5:max5
        LandauM = m-min5;
        % <phi_5_m|H|phi_9n>
        for n = min9:max9
            LandauN = n-min9;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma5-gamma2)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
        % <phi_5_m|H|phi_11n>
        for n = min11:max11
            LandauN = n-min11;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma2+gamma5)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
    end
    
    
    %<phi_6_m|H|phi> is zero except <phi_6_m|H|phi_10n> and <phi_6_m|H|phi_12n>
    
    for m = min6:max6
        LandauM = m-min6;
        % <phi_6_m|H|phi_9n>
        for n = min10:max10
            LandauN = n-min10;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma2 - gamma5)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
        % <phi_6_m|H|phi_12n>
        for n = min12:max12
            LandauN = n-min12;
            Ht(m,n) = Ht(m,n) + 1/2*(gamma5 + gamma2)*delta(LandauM,LandauN);
            Ht(n,m) = Ht(m,n);
        end
        
    end


%**********************************************************
    
% load('0223test.mat');
% 
% diffH=Ht-Htt;
% maxelem=max(max(diffH))

    [Q,D] = eig(Ht);
%     clear Ht;
    E(step,:)=[B;diag(D)];
 end

% filename = strcat('Kp/Kp_EIGEN_G4_LL-E_',num2str(Eadd),'-N_',int2str(N),'-g_0=',num2str(gamma0),'-g_1=',num2str(gamma1),'-g_2=',num2str(gamma2),'-g_3=',num2str(gamma3),'-g_4=',num2str(gamma4),'-g_5=',num2str(gamma5),'-DAB=',num2str(deltaAB),'.mat')
%  save(filename);

figure1=figure;

axis([Bmin Bmax Emin Emax]);
hold on

for i=2:5+4*N
    plot(E(:,1),E(:,i),'r','linewidth',0.5);
end

for i= 6+4*N:9+4*N
    plot(E(:,1),E(:,i),'k','linewidth',0.5);
end

for i= 10+4*N:13+8*N
    plot(E(:,1),E(:,i),'b','linewidth',0.5);
end


hold off

xlabel('B(T)')
ylabel('E(meV)')

str = strcat('r1=',num2str(gamma1),'r3=',num2str(gamma3),'r4=',num2str(gamma4),'r2=',num2str(gamma2),',r5=', num2str(gamma5),...
                 ',d=',num2str(deltaAB),',d_middle=', num2str(delta_middle),'E=',num2str(Eadd), ',N=', num2str(N),',step=', num2str(stepsB));

%str = strcat('N=',int2str(N));
            
annotation('textbox','String',str,'FitBoxToText','on');
% figname = strcat('r1BLG=',num2str(r1BLG),',r1intr=',num2str(r1intr),',r3=',num2str(gamma3),',r4=',num2str(gamma4),',r2=',num2str(gamma2),',r5=', num2str(gamma5),...
%                    ',d=',num2str(deltaAB),',d_middle=', num2str(delta_middle),'E=',num2str(Eadd),',N=', num2str(N),',step=', num2str(stepsB),'.jpg');
% saveas(figure1,figname);





