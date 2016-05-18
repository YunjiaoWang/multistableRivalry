epsilon = 0.01; % slow parameter
params(1) = epsilon;
alpha = 0.5; % Maximum excitation rate between cells on the same row
params(2) = alpha ;
beta = 0.5; % Maximum excitation rate between cells on the same diagonal:
params(3) = beta ;

%delta = 1.5; % 1, Inhibition rate between cells on the same column
delta = 1;
params(4)=delta;
%gamma = 1; % 1, Fatigue rate
gamma=0.9;
%gamma = 0; % Fatigue rate
params(5) = gamma;
I = 1; % input
params(6) = I;
%%% parameter for upper level
delta1 = 0.3;
 params(7) = delta1;
delta2 = 0.3;
params(8) = delta2 ;
gamma1 = 0.9;
params(9) = gamma1;
epsilonS = epsilon;
params(10)=epsilonS;
phi=1;
params(11) = phi;
m = 0.1;
params(12) = m;

