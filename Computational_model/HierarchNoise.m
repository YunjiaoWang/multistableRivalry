function dy=HierarchNoise( y, noise,params)

%% Model Description.%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% network model at the first level%%%%
% 1 <-->2
% -     -
% |     |
% -     -
% 3 <-->4
% Excitation between nodes on the same row: nodes 1 and 2, nodes 3 and 4.
% Excitation between nodes on the same diagonal: nodes 1 and 4, nodes 2 and 3
% Inhibition between  nodes on the same column: nodes 1 and 3, nodes 2 and 4 

%% % network model at the second level %%%
% Four nodes representing populations of four percepts (P1, P2, P3 and P4) inhibite each other.
% Signal from nodes 1 and 2 excites P1, nodes 3 and 4 excites P2, nodes 1
% and 4 excites P3, and nodes 2 and 3 excites P4.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = params(1);
alpha = params(2);
beta = params(3);
delta = params(4);
gamma = params(5);
I = params(6);
%% % parameter for upper level
delta1 = params(7);
delta2 = params(8);
gamma1 = params(9);
epsilonS = params(10);
%%
dy = zeros(16, 1);
xe = zeros(4, 1); % activity variables at the first level
xh = zeros(4, 1); % fatigue variable at the first level
 E = zeros(4, 1); % activity variables at the 2nd level
 H = zeros(4, 1); % fatigue variable at the 2nd level

%%   
for k = 1:4 
  xe(k) = y(2*k - 1);
  xh(k) = y(2*k );
   E(k) = y(8 + 2*k - 1);
   H(k) = y(8+ 2*k);
end


%% %%%% Downstream %%%%%%%%%%%%%%%%%%
% equations of x1= (xe(1), xh(1))
dy(1) = -xe(1) + G(alpha*xe(2) + beta * xe(4) - delta*xe(3) - gamma*xh(1) + I+noise(1)); 
dy(2)= epsilon*(- xh(1) + xe(1));

% equations of x2= (xe(2), xh(2))
dy(3) = -xe(2) + G(alpha*xe(1) + beta  * xe(3) - delta*xe(4) - gamma*xh(2) + I +noise(2)); 
dy(4)= epsilon*(- xh(2) + xe(2));

% equations of x3= (xe(3), xh(3))
dy(5) = -xe(3) + G(alpha*xe(4) + beta  * xe(2) - delta*xe(1) - gamma*xh(3) + I+noise(3)); 
dy(6)= epsilon*(- xh(3) + xe(3));


% equations of x4 = (xe(4), xh(4))
dy(7) = -xe(4) + G(alpha*xe(3) + beta * xe(1) - delta*xe(2) - gamma*xh(4) + I  +noise(4)); 
dy(8)= epsilon*(- xh(4) + xe(4));

%% %%%%%%%%%% upper Stream %%%%%%%%%%%%%%%
% % four percepts
% percept 1 associated to the 12 patches, equations of (E1, H1)
 dy(9) = -E(1) + G(xe(1)*xe(2) - delta1*E(2) - delta2*E(3) - delta2*E(4) - gamma1*H(1)+noise(9));
 dy(10) = epsilonS*(E(1) - H(1));
%  percept 2 associated to the 34 patches, equations of (E2, H2)
 dy(11) = -E(2) + G((xe(3)*xe(4)) - delta1*E(1) - delta2*E(3) - delta2*E(4) - gamma1*H(2)+noise(10));
dy(12) = epsilonS*(E(2) - H(2));
 
% %  percept 3 associated to the 14 patches, equations of (E3, H3)
 dy(13) = -E(3) + G((xe(1)*xe(4)) - delta1*E(4) - delta2*E(1) - delta2*E(2) - gamma1*H(3)+noise(11));
dy(14) = epsilonS*(E(3) - H(3));

% % percept 4 associated to the 23 patches, equations of (E4, H4)
 dy(15) = -E(4) + G((xe(2)*xe(3)) - delta1*E(3) - delta2*E(1) - delta2*E(2) - gamma1*H(4)+noise(12));
dy(16) = epsilonS*(E(4) - H(4));


%%%%%%%%%%%%%%
function y = G(x)
%% gain function
% Gain funtion usually takes sigmoidal form. The specific function can either be smooth or nonsmooth
% Here, we take the smooth form.

%%%%% Smooth function%%%%%%%
a = 0.4;
b = 0.2;
c = 10;
y = 2*a/(1+exp(-c*(x-b)));
