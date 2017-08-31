function f = phi_fn(z,d,Phi_L,Phi_R,eV)
%%--- Parameters ---%%
Phi_I = 0; % electron-affinity of insulator [eV]

%%--- Function ---%% %[eV]
f_simple = (Phi_L - Phi_I) + ((Phi_R-Phi_L) - eV).*z./d;
f_simple(1) = 0; f_simple(end) = -eV; % border conditions
f = f_simple;
