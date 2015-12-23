
%% parameters for the 1st level network model
epsilon = 0.01; % slow parameter
params(1) = epsilon;
alpha = 0.5; % Maximum excitation rate between cells on the same row
params(2) = alpha ;
beta = 0.5; % Maximum excitation rate between cells on the same diagonal:
params(3) = beta ;

delta = 0.9; % 1, Inhibition rate between cells on the same column
params(4)=delta;
gamma = 0.9; % 0.9, Fatigue rate
params(5) = gamma;
I = 0.7; % input
params(6) = I;


%% parameters for 2nd level network model
delta1 = 0.5;
 params(7) = delta1;
delta2 = 0.5;
params(8) = delta2 ;
gamma1 = 1;
params(9) = gamma1;
epsilonS = epsilon;
params(10)=epsilonS;
