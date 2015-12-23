function noise = ColorNoise(Dt,  N, m, sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Ornstein-Uhlenbeck process to produce colored noise:
% The noise n(t) is an OU process with zeros means and 
% deviation sigma. i.e. n(t) is a solution to
%   
% n'(t) = -\frac{n}{taus} + sigma sqrt{\frac{2}{taus}} xi(t)
% 
% where xi(t)is a white noise process with zero mean and variance
% t-t' for time t'<t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (See the equation (A1) in the Appendix of Moreno et al. 07 paper)
% Also I follows the algorithm that was proposed in 
% 'An Algorithmic: Introduction to
% Numerical Simulation of Stochastic Different ial Equations' 
%( by Desmond J. Higham,  SIAM REVIEW, 2001 
% Vol. 43,No . 3,pp . 525-546) for generating the noise sequence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 2;
dt = Dt/R;
M = R*N; 
dw = sqrt(dt)*randn(M, m);
tau = 1;
taus = 20; % milisecond
a = 1/taus;
%sigma = 0.03;
Xzero = sigma*randn(1,m); 
%Xzero = zeros(1, m);
noise = zeros(N, m);

Xtemp = Xzero;
for j= 2:N
    Winc = sum(dw(R*(j-1)+1:R*j, :));
    Xtemp = Xtemp -tau* Dt*a*Xtemp +tau* sigma*sqrt(2*a)*Winc;
    noise(j,:) = Xtemp;
end

