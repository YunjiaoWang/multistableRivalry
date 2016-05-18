%% seaching for the dominance time for each percept
function [Domp] = FindDominancePeriods(T, Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domp:1st col - percept #, 2nd col - duration 
% other two colums are starting and ending time.
% i.e. 
% Domp(:,1) -- percept
% Domp(:, 2) -- period
% Domp(:, 3)-- starting time
% Domp(:, 4)-- ending time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1 = Y(:, 9);
P2 = Y(:, 11); 
P3 = Y(:, 13);
P4 = Y(:, 15);
L = length(T); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13 = percept1; 24 = percept2;  14 = percept3; 23 = percept4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er = 0.05;
Domp = zeros(100,4);
tStart = 0 ;
tEnd = 0;

i = 1 ;
p = 0 ;
% Record the beginnings and endings of each dominanting period
while i<=L 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Searching for percept1  dominating
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if (P1(i) > P2(i)+er) && (P1(i) > P3(i) + er)&& (P1(i) > P4(i) + er)
        %u = u+1; 
        tStart = T(i);
        % searching for the end
        i = i+1;
        while i<=L && (P1(i) > P2(i)+er) && (P1(i) > P3(i) + er)&& (P1(i) > P4(i) + er)
            i = i+1;
        end
        tEnd = T(i-1);
        p = p+1; 
        Domp(p, 1) = 1;
        Domp(p, 2) = tEnd - tStart;
        Domp(p, 3) = tStart;
        Domp(p, 4) = tEnd;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % searching for percept2(24) dominating
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (P2(i) > P3(i)+er) && (P2(i) > P1(i) + er)&& (P2(i) > P4(i) + er)
        tStart = T(i);
        % searching for the end
        i = i+1;
        while i<=L && (P2(i) > P3(i)+er) && (P2(i) > P1(i) + er)&& (P2(i) > P4(i) + er)
            i = i+1;
        end
        tEnd = T(i-1);
        p = p+1; 
        Domp(p, 1) = 2;
        Domp(p, 2) = tEnd - tStart; 
        Domp(p, 3) = tStart;
        Domp(p, 4) = tEnd;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % searching for percept3(14) dominating
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif  (P3(i) > P2(i)+er) && (P3(i) > P1(i) + er)&& (P3(i) > P4(i) + er)
        tStart = T(i);
        % searching for the end
        i = i+1;
        while i<=L && (P3(i) > P2(i)+er) && (P3(i) > P1(i) + er)&& (P3(i) > P4(i) + er)
            i = i+1;
        end
        tEnd = T(i-1);
        p = p+1; 
        Domp(p, 1) = 3;
        Domp(p, 2) = tEnd - tStart; 
        Domp(p, 3) = tStart;
        Domp(p, 4) = tEnd;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % searching for percept4(23) dominating
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    elseif  (P4(i) > P2(i)+er) && (P4(i) > P1(i) + er)&& (P4(i) > P3(i) + er)
        tStart = T(i);
        % searching for the end
        i = i+1;
        while i<=L && (P4(i) > P2(i)+er) && (P4(i) > P1(i) + er)&& (P4(i) > P3(i) + er)
            i = i+1;
        end
        tEnd = T(i-1);
        p = p+1; 
        Domp(p, 1) = 4;
        Domp(p, 2) = tEnd - tStart; 
        Domp(p, 3) = tStart;
        Domp(p, 4) = tEnd;
    else
        i=i+1;
    end
end
Domp = Domp(1:p,:);
