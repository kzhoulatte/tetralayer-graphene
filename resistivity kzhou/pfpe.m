function  y = pfpe(E,Ef,KT)

n = exp((E-Ef)/KT);
n_1 = exp(-(E-Ef)/KT);

y = 1/KT * 1/(n+n_1+2);
end