%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculation of the Quantum Transmission Probability %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Javier Carrasco / javier.carrasco@ug.uchile.cl
% -------------------------------------------------------
% DESCRIPTION:
%
% This program calculates the transmission probability, T(E), as a
% function of energy, E, for many different values of bias voltage
% applied between the contacts.
%
% The intention is to calculate T(E) for MIM junctions,
% so the potential energy at the extremes represent the conduction band
% energy levels of the respective metals and the potential energy inbetween
% corresponds to the conduction band energy level of the insulator (which
% might change a lot depending on the configuration chosen), so I treat
% this region as an arbitrary energy potential shape region (it could also
% be multiple insulators, for example).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%--- Parameters ---%%
c = 2.99792458e8; % speed of light [m/s^2]
h_bar = 6.58211899e-16; % reduced Planck's constant [eV*s]
m_e = 0.5109989461*10^6/(c^2); % rest electron mass [eV*s^2/m]
Phi_L = 4; % left-electrode's electron-affinity [eV]
Phi_R = 4; % right-electrode's electron-affinity [eV]
d = 0.7e-9; % insulator(s) width [m]
epsilon_0 = 8.8541878176e-12; % permitivity of free space [C^2/(N*m^2)]
                              %                           (SI)
epsilon_r = 1; % relative permitivity of the insulator [ ]
e_charge = 1.6021766208e-19; % electron charge [C]

d_c = h_bar/((2*m_e*Phi_L)^(1/2)); % scaling factor [m]?? review units

for V = -2:1:2 % bias voltage [V]

eV = V; % [eV] (if I multiply by e, then it is in [Joules], not [eV])

dE = 0.01;
z_relative_step = 0.01;
dz = z_relative_step*d;

T_List = [];

    for E = 0:dE:15

        %%---Functions---%%
        z = 0:dz:d; % n+1 boundaries, such that n*0.01*d=d => n=10^2
        % [m]
        z = z/d_c; % scale to d_c

        slicesNumber = length(z); % n+1=10^2+1 in this case
        MatrixSize = slicesNumber*2;
        n = slicesNumber-1; %n; for easier manipulation later
        phi = phi_fn(z,d/d_c,Phi_L,Phi_R,eV); % Potential Energy Barrier; [eV]

        k_0 = sqrt(2*m_e*E)/h_bar; % k_0 def.
        % [m^(-1)]
        % These (the later) are the units for all the k's... I'll change this later
        k_nPlus1 = sqrt(2*m_e*(E+eV))/h_bar; % k_(n+1) def.

        k = (2*m_e.*(phi-E)).^(1/2)/(h_bar); % n+1 k-values, but only internal values
        % and the k_(n+1) is wrong (since it is not internal
        %[i.e. not inside the barrier])
        k(n+1) = k_nPlus1; %Correct the k_(n+1) value
        k = [0,k]; %Add a 0th value (l=0)
        k(1) = k_0; %Correct the 0th value
        % k-vector with all the k_l, l=0,...,n+1 values; fully written
        % Remeber: indeces in Matlab start from 1, not 0

        %%%% k=(10^9)*k; % this redefines the units as [nm^(-1)]

        k = k*(d_c); % scale to k_c = (1/d_c)

        % Note: lenght(k)=n+2

        %%---Matrices---%%
        X=zeros(MatrixSize,1); % column-vector of length 2*(n+1); with zeros
        X([1,2])=[-1;-1i*k_0*d_c]; % redefine first 2 values;
        %here X is completely defined;
        %units: [1];[nm^(-1)];[1];......;[nm^(-1)]  (change in each row)

        M=zeros(MatrixSize); % 2*(n+1) by 2*(n+1) matrix of zeros

        M(1,1)=1; M(2,1)=-1i*k_0*d_c; % M_0^(-) written

        M(MatrixSize-1,MatrixSize)=-exp(1i*k(n+1+1)*d/d_c);
        M(MatrixSize,MatrixSize)=-1i*(k(n+1+1))*exp(1i*k(n+1+1)*d/d_c);
        % -M_n^(+) written

        s=1; %RowNumber_index
        l=0; %MatrixBlock_index
        while l<=n-1
            M(s,s+1)=-exp(k(l+1+1)*z(l+1));
            M(s,s+2)=-exp(-k(l+1+1)*z(l+1));
            M(s+1,s+1)=-(k(l+1+1))*exp(k(l+1+1)*z(l+1));
            M(s+1,s+2)=(k(l+1+1))*exp(-k(l+1+1)*z(l+1));
            %-M_(l+1)^(+) written
            l=l+1;
            %l-index updated
            M(s+2,s+1)=exp(k(l+1)*z(l+1));
            M(s+2,s+2)=exp(-k(l+1)*z(l+1));
            M(s+3,s+1)=(k(l+1))*exp(k(l+1)*z(l+1));
            M(s+3,s+2)=-(k(l+1))*exp(-k(l+1)*z(l+1));
            %M_l^(-) written
            s=s+2; %update s-index
        end

        %Units of M: if s is even => that row has [nm^(-1)];
        %            if s is odd => that row has [1];

        % Therefore, as M*A=X, the (yet to be defined, but it is
        % the variable vector to solve... the coefficients) A vector is
        % dimensionless, i.e. it has [1]

        %%--- Solve Matrix System ---%%
        A = (M^(-1))*X;

        %To calculate the transmission probability, the only coeff. we are
        %interested in is the A_(n+1)^(+)... which is the last one of the A vector
        %Remember, in the theory we have set A_0^(+)=1

        A_final = A(end); %gets the last coordinate-value of the A vector

        %%--- Obtain the Transmission Probability ---%%
        T = (k_nPlus1/k_0)*(abs(A_final)^2);

        %T; %print T to check values in the console
        T_List = [T_List T];

    end

    E=0:dE:15;

    figure(1)
    plot(E,T_List)
    xlabel('E [eV]'); ylabel('T(E)')
    title('Transmission Probability v/s Energy')
    grid on
    hold on

    %%% To plot the phi vs z/d_c:
    %hold off
    figure(2)
    plot(z,phi)
    xlabel('z [d_c]'); ylabel('\phi(z) [eV]')
    title('Potential Barrier v/s Position')
    %axis([0 d/d_c -6 5])
    hold on
    grid on

end

%Prepare legends for figures:
V = -2:1:2;
Length_V = length(V);

Legend_cell = cell(Length_V,1);
for iter = 1:Length_V
    Legend_cell{iter} = strcat('V = ', num2str(V(iter)));
end
%-----------------------------

%Add legends to figures:
figure(1)
legend(Legend_cell,'Location','northeastoutside')
figure(2)
legend(Legend_cell,'Location','northeastoutside')
%------------------------

%Save figures:
figure(1)
print(strcat('T_step',num2str(z_relative_step),'_varV.eps'),'-depsc')
figure(2)
print(strcat('phi_step',num2str(z_relative_step),'_varV.eps'),'-depsc')
%-------------

