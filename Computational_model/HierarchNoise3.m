function dy=HierarchNoise3( y, noise,params)
%% Model Description.
% network model 
% 1 <-->2
% -     -
% |     |
% -     -
% 3 <-->4
% Excitation between cells on the same row: cells 1 and 2, cells 3 and 4.
% Excitation between cells on the same diagonal: cells 1 and 4, cells 2 and 3
% Inhibition between cells on the same column: cells 1 and 3, cells 2 and 4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%epsilon = 0.01; % slow parameter
epsilon = params(1);
%alpha = 0.48; % Maximum exicitation rate between cells on the same row
alpha = params(2);
%beta = 0.5; % Maximum exicitation rate between cells on the same diagonal:
beta = params(3);

%delta = 1; % Inhibation rate between cells on the same column
delta = params(4);
%gamma = 1; % Fatigue rate
gamma = params(5);
%I = 0.51; % input
I = params(6);
%%% parameter for upper level
%delta1 = 0.5;
delta1 = params(7);
%delta2 = 0.4;
delta2 = params(8);
%gamma1 = 1;
gamma1 = params(9);
%epsilonS = epsilon;
epsilonS = params(10);
phi = params(11);
%m = params(12);
%%
dy = zeros(16, 1);
xe = zeros(4, 1); % activity variable
xh = zeros(4, 1); % fatigue
 E = zeros(4, 1); % adaptation of excitation coupling between row cells 
 H = zeros(4, 1); % adaptation of excitation coupling between diagonal cells

%%   
for k = 1:4 
  xe(k) = y(2*k - 1);
  xh(k) = y(2*k );
   E(k) = y(8 + 2*k - 1);
   H(k) = y(8+ 2*k);
end

%% Downstream
% equations of x1= (xe(1), xh(1))
%dy(1) = -xe(1) + G(alpha*xe(3) + beta * xe(4) - delta*xe(2) - gamma*xh(1) + I + sigma*(E(1)+E(3))*xe(3)); 
% % Without Noise
% dy(1) = -xe(1) + G(alpha*xe(2) + beta * xe(4) - delta*xe(3) - gamma*xh(1) + I); 
% With noise
dy(1) = -xe(1) + G(alpha*xe(2) + beta *xe(4) - delta*xe(3) - gamma*xh(1) + I+noise(1)); 
%dy(2)= epsilon*(- xh(1) + xe(1)+noise(5));
dy(2)= epsilon*(- xh(1) + xe(1));

% equations of x2= (xe(2), xh(2))
% %Without noise
% dy(3) = -xe(2) + G(alpha*xe(1) + beta  * xe(3) - delta*xe(4) - gamma*xh(2) + I ); 
%With noise
dy(3) = -xe(2) + G(alpha*xe(1) + beta  * xe(3) - delta*xe(4) - gamma*xh(2) + I +noise(2)); 
% dy(3) = -xe(2) + G(alpha*xe(4) + beta  * xe(3) - delta*xe(1) - gamma*xh(2) + I + sigma*(E(2)+E(4))*xe(4)); 
%dy(4)= epsilon*(- xh(2) + xe(2)+noise(6));
dy(4)= epsilon*(- xh(2) + xe(2));

% equations of x3= (xe(3), xh(3))
% %Without noise
% dy(5) = -xe(3) + G(alpha*xe(4) + beta  * xe(2) - delta*xe(1) - gamma*xh(3) + I); 

%With noise
dy(5) = -xe(3) + G(alpha*xe(4) + beta * xe(2) - delta*xe(1) - gamma*xh(3) + I+noise(3)); 
%dy(5) = -xe(3) + G(alpha*xe(1) + beta  * xe(2) - delta*xe(4) - gamma*xh(3) + I + sigma*(E(1)+E(4))); 
%dy(6)= epsilon*(- xh(3) + xe(3)+noise(7));
dy(6)= epsilon*(- xh(3) + xe(3));


% equations of x4 = (xe(4), xh(4))
% % Without noise
% dy(7) = -xe(4) + G(alpha*xe(3) + beta * xe(1) - delta*xe(2) - gamma*xh(4) + I  +noise(4)); 

% With noise
dy(7) = -xe(4) + G(alpha*xe(3) + beta* xe(1) - delta*xe(2) - gamma*xh(4) + I  +noise(4)); 
%dy(7) = -xe(4) + G(alpha*xe(2) + beta * xe(1) - delta*xe(3) - gamma*xh(4) + I + sigma*(E(2)+E(3))); 
%dy(8)= epsilon*(- xh(4) + xe(4)+noise(8));
dy(8)= epsilon*(- xh(4) + xe(4));

% %%%%%%%%%% upper Stream %%%%%%%%%%%%%%%
% % four percepts
% % percept 1 associated to the 12 patches, equations of (E1, H1)
% with no noise 
%dy(9) = -E(1) + G(xe(1)*xe(2) - delta1*E(2) - delta2*E(3) - delta2*E(4) - gamma1*H(1));
% %% with Noise
 dy(9) = -E(1) + G(phi*(xe(1)*xe(2)) - delta1*E(2) - delta2*E(3) - delta2*E(4) - gamma1*H(1)+noise(5));
%  %dy(9) = -E(1) + G((xe(1)+xe(2)) - delta1*E(2) - delta2*E(3) - delta2*E(4) - gamma1*H(1));
% With no noise
 dy(10) = epsilonS*(E(1) - H(1));
%   % With noise
 %dy(10) = epsilonS*(E(1) - H(1)+noise(13));

% % percept 2 associated to the 34 patches, equations of (E2, H2)
% without noise 
%dy(11) = -E(2) + G((xe(3)*xe(4)) - delta1*E(1) - delta2*E(3) - delta2*E(4) - gamma1*H(2));
% % with noise 
 dy(11) = -E(2) + G(phi*(xe(3)*xe(4)) - delta1*E(1) - delta2*E(3) - delta2*E(4) - gamma1*H(2)+noise(6));
% Without noise 
dy(12) = epsilonS*(E(2) - H(2));
% % % With noise 
 %dy(12) = epsilonS*(E(2) - H(2) + noise(14));
% % 
% %  percept 3 associated to the 14 patches, equations of (E3, H3)
% Without noise
%dy(13) = -E(3) + G((xe(1)*xe(4)) - delta1*E(4) - delta2*E(1) - delta2*E(2) - gamma1*H(3));
%% With noise
 dy(13) = -E(3) + G(phi*(xe(1)*xe(4)) - delta1*E(4) - delta2*E(1) - delta2*E(2) - gamma1*H(3)+noise(7));
%   %Without noise
dy(14) = epsilonS*(E(3) - H(3));
%   %With noise
 %dy(14) = epsilonS*(E(3) - H(3) + noise(15));

% % percept 4 associated to the 23 patches, equations of (E4, H4)
% Without noise 
%dy(15) = -E(4) + G((xe(2)*xe(3)) - delta1*E(3) - delta2*E(1) - delta2*E(2) - gamma1*H(4));
% %% With noise 
 dy(15) = -E(4) + G(phi*(xe(2)*xe(3)) - delta1*E(3) - delta2*E(1) - delta2*E(2) - gamma1*H(4)+noise(8));
 % Without noise 
dy(16) = epsilonS*(E(4) - H(4));
% % With noise 
  %dy(16) = epsilonS*(E(4) - H(4) + noise(16));

%%%%%%%%%%%%%%
function y = G(x)
%% gain function
% Gain funtion usually takes sigmoidal form. The specific function can either be smooth or nonsmooth
% Here, we take the smooth form.

%%%%% Smooth function%%%%%%%
a = 0.4;
b = 0.5;
c = 10;
y = 2*a/(1+exp(-c*(x-b)));
%%%%%%% Nonsmooth function %%%%%%%%
%y = (x>0.5); 

% function y = f(x)
% y = (abs(x) - 1>0);
