%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This file is to create the four main numerical results in the paper:
% Color saturation enhances interocular grouping in perceptual multistability
% By Yunjiao Wang, Alain Jacot-Guillarmod, Claudia Pedroza, Haluk Ogmen, Kresimir Josic, Zachary Kilpatrick
% 2016.
% More specifically, the file will generate the following graphs:
% 1. Predominance of grouped percepts vs grouping strength beta
% 2. Dominance durations of grouped and single-eye percepts vs beta
% 3. The ratios of number of visits to grouped percepts  vs beta
% 4. The transition probabilities  vs beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


parameters1;
rand('state', 200);
Dt = 0.1; %ms
tstart = 0;
tend = 2000;
tspan = tstart:Dt:tend;
T = tspan;
N = length(tspan);

nSeq = 10; % number of time series under each given beta value.

b0 = 0.2;
b1 = 0.8;
db = 0.01;

nb = floor((b1-b0)/db);

% set alpha value
params(2) = 0.8; 


%%%%%%%%%%%%%%%%%%%%%
%% Simulating data.%%
%%%%%%%%%%%%%%%%%%%%%
Dom_Single_mean = zeros(nb, 1); % average dominance duration of single eye percepts
Dom_Grouped_mean = zeros(nb, 1); % average dominance duration of single eye percepts
nVisit_Group_mean = zeros(nb, 1); % the ratio of the number of visits to grouped percepts  n_grouped /(n_grouped + n_single)
Dom_ratiosGroup_mean = zeros(nb, 1); % Predominance of grouped percepts: T_Grouped/ (T_Grouped + T_single)

Dom_Single_std = zeros(nb, 1); % standard deviation of dominance duration of single eye percepts
Dom_Grouped_std = zeros(nb, 1); % standard deviation of dominance duration of single eye percepts
nVisit_Group_std = zeros(nb, 1); % standard deviation of the ratio of the number of visits to grouped percepts
Dom_ratiosGroup_std = zeros(nb, 1); %standard deviation of predominance: T_Grouped/ (T_Grouped + T_single)

%%  transition
Single_mean = zeros(nb, 1); % mean value of transition probability from single to single
Grouped_mean = zeros(nb, 1); % mean value of transition probability from grouped to grouped
Single_to_grouped_mean = zeros(nb, 1); % mean value of transition probability from single to grouped
grouped_to_single_mean = zeros(nb, 1); % mean value of transition probability from grouped to single


Single_std = zeros(nb, 1); % standard deviation of transition probability from single to single
Grouped_std = zeros(nb, 1); % standard deviation of transition probability from grouped to grouped
Single_to_grouped_std = zeros(nb, 1); % standard deviation of transition probability from single to grouped
grouped_to_single_std = zeros(nb, 1); % standard deviation of transition probability from grouped to single


% params(2)=0.5;
betaSeq = b0:db:b1;
cut = 4; % the theshold that cut off the piece with short period.
%%
for b=1:nb+1  
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    params(3) = betaSeq(b);    
    
    
    %%%% For computing average dominance durations of single-eye/grouped perceptss %%%%%
    dompSingle = zeros(nSeq*50, 2); % initialize
    nSingle = 0;
    dompGroup = zeros(nSeq*50, 2);% initialize
    nGroup = 0;
    
    %%%%% for computing ratios%%%%%%%%%
    % for each trial, we need to find the ratios of the number of visits
    % to grouped percepts, and ratios of the total dominance of grouped
    % percepts.
    nVisRatioGrouped = zeros(nSeq,1); %the ratios of the number of visits to grouped percepts
    TotDomRatioGrouped = zeros(nSeq,1); %ratios of the total dominance of grouped % percepts.
    
    
    %%%%%%%% for computing transition probability %%%%%%%%%%%%%%%%%%%%  
     trans_ss = zeros(nSeq,1);
     trans_sg = zeros(nSeq,1);
     trans_gg = zeros(nSeq,1);
     trans_gs = zeros(nSeq,1);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate a sequence of time series and compute each aimed quantity.
     for u = 1:nSeq
        %%%%% Add colored noise;
        % m  -- number of series noise generated.
        m = 16;
        sig = 0.07; %%% standard deviation.
        noise = ColorNoise(Dt, N, m, sig);
        %%%%%%%%%%%%%%%%%%%
        %%% Initialize.
        yzero = rand(16, 1);
        %%%% Euler method
        Y = zeros(N, 16);
        Y(1,:) = yzero';
        for i  = 2: N
            yn = Y(i-1, :)';
            fn = HierarchNoise(yn, noise(i, :),params);
            ynplus1 = yn + Dt*fn;
            Y(i, :) = ynplus1';
        end
        %Analysis
        %%% throw away the first second data %%%
        Y= Y(100:end,:); %%
        T1 = T(100:end);
        
        %%% throw away the dominance duration that is less than the cut off
        %%% value
        domp = FindDominancePeriods(T1, Y);        
        domDur = domp(:, 2);
        [r,c] = find(domDur>cut);
        domp = domp(r, :);
        
        %%% computing transition ratio of the resulting time series.
        len = length(domp);
        A = transMatrix(domp);
        
        A11 = A(1:2, 1:2);
        A22 = A(3:4, 3:4);
        A21 = A(3:4, 1:2);
        A12 = A(1:2, 3:4);
        
        trans_ss(u) = sum(A11(:)) / (sum(A11(:)+sum(A12(:))));
        trans_sg(u) = 1-  trans_ss(u);
        trans_gg(u) = sum(A22(:)) / (sum(A21(:)+sum(A22(:))));
        trans_gs(u) = 1 - trans_gg(u) ;
        
        %%% find rows of single-eype percepts %%%
        perc = domp(:,1);
        [r,c]=find(perc<3);
        singleDom = domp(r,1:2);
         %%% find rows of grouped percepts %%%
        [r,c]=find(perc>2);
        groupedDom = domp(r,1:2);
        
        %%% compute the ratios of the number of visits to grouped percepts 
        numVistGrouped = length(groupedDom(:,1));
        numVistSingle = length(singleDom(:,1));
        nVisRatioGrouped(u) =  numVistGrouped/(numVistGrouped+numVistSingle);
        
        %%%% compute predominance of grouped percepts %%%
        DomGrouped = sum(groupedDom(:,2));
        DomSingle = sum(singleDom(:,2));
        TotDomRatioGrouped(u) = DomGrouped/(DomGrouped+DomSingle);
        
        %%%% put all of the dominance durations of single-eye perc. in one vector
        if (numVistSingle>0)
            dompSingle(nSingle+1:nSingle+numVistSingle,:) = singleDom;
            nSingle = nSingle+numVistSingle;
        end
        %%%% put all of the dominance durations of grouped perc. in one vector
        if(numVistGrouped>0)
            dompGroup(nGroup+1:nGroup+numVistGrouped,:) = groupedDom;
            nGroup = nGroup + numVistGrouped;
        end
     end % End of the nSeq.
     
     %% Find mean and std of each quantity.
     
     %%%% throw away empty elements %%%
     dompSingle=dompSingle(1:nSingle,:);
     dompGroup = dompGroup(1:nGroup,:);
     
     % compute mean and variance of single-eye percepts
     Dom_Single_mean(b) = mean(dompSingle(:,2)); 
     Dom_Single_std(b) = std(dompSingle(:,2)); 
     
     % compute mean and variance of grouped percepts
     Dom_Grouped_mean(b) = mean(dompGroup(:,2)); 
     Dom_Grouped_std(b) = std(dompGroup(:,2)); 
     
     % compute mean and std of the number of visits to grouped percepts
     nVisit_Group_mean(b) = mean(nVisRatioGrouped);
     nVisit_Group_std(b) = std(nVisRatioGrouped);
     
      % compute mean and std of the predominance of grouped percepts
     Dom_ratiosGroup_mean(b) = mean(TotDomRatioGrouped); %
     Dom_ratiosGroup_std(b) = std(TotDomRatioGrouped); %
     
     % compute meand and std of transition probability
     Single_mean(b) = mean(trans_ss); % single to single
     Single_std(b) = std(trans_ss); % single to single
     
     Grouped_mean(b) = mean(trans_gg); % grouped to grouped
     Grouped_std(b) = std(trans_gg); % grouped to grouped
     
     Single_to_grouped_mean(b) = mean(trans_sg); % single to grouped
     Single_to_grouped_std(b) = std(trans_sg); % single to grouped
     
     grouped_to_single_mean(b) = mean(trans_gs); % grouped to single
     grouped_to_single_std(b) = std(trans_gs); % grouped to single
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of beta loop


%% Plot the Analysis

B = b0:db:b1;
B = B(1:nb+1);

%%% plot predominance of grouped percepts %%%
figure(1)
y1 = Dom_ratiosGroup_mean - Dom_ratiosGroup_std;
y2 = Dom_ratiosGroup_mean + Dom_ratiosGroup_std;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.1 0.1 1]);
hold on
plot(B, Dom_ratiosGroup_mean, 'r', 'LineWidth', 3);
ylabel('Predominance', 'FontSize',30);
title('Grouped percepts')
set(gca, 'FontSize',24)
hold off

%%% plot the ratios of the number of visits to grouped percepts %%%
figure(2)
y1 =  nVisit_Group_mean -  nVisit_Group_std;
y2 =  nVisit_Group_mean +  nVisit_Group_std;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 1]);
hold on
plot(B, nVisit_Group_mean, 'k', 'LineWidth', 3)
%title('Average Ratios of visits to Grouped percepts', 'FontSize',20)
xlabel('\beta', 'FontSize',30);
ylabel('Ratio', 'FontSize',30);
title('Num. of visits to grouped perc.')
set(gca, 'FontSize',24)
hold off



%%% plot the mean and std of dominance durations of grouped percepts/single-eye percepts %%%
figure(3)
y1 = (Dom_Grouped_mean - Dom_Grouped_std)/10;
y2 = (Dom_Grouped_mean + Dom_Grouped_std)/10;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 1]);
hold on
y1 = (Dom_Single_mean - Dom_Single_std)/10;
y2 = (Dom_Single_mean + Dom_Single_std)/10;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 1]);
p=plot(B, Dom_Grouped_mean/10, B, Dom_Single_mean/10,'LineWidth', 3);
legend(p,'grouped to grouped', 'single to single')
ylabel('Dom. Duration(sec)', 'FontSize',30);
xlabel('\beta', 'FontSize', 30)
set(gca, 'FontSize',24)
hold off

%%% plot the mean and std of transition probabilities %%%
figure(4)

[r,c] = find(~isnan(Grouped_mean));
y1 = Grouped_mean(r) - Grouped_std(r);
y2 = Grouped_mean(r) + Grouped_std(r);
B1 = B(r);
fill_between(B(r), y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.98 0.97 0.97]);
hold on

[r,c] = find(~isnan(Grouped_mean));
y1 = Single_mean - Single_std;
y2 = Single_mean + Single_std;
y1 = y1(r);
y2 =y2(r);
fill_between(B(r), y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]);

y1 = Single_to_grouped_mean - Single_to_grouped_std;
y2 = Single_to_grouped_mean + Single_to_grouped_std;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 1]);

y1 = grouped_to_single_mean - grouped_to_single_std;
y2 = grouped_to_single_mean + grouped_to_single_std;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.1 0.1 1]);

p=plot(B,Grouped_mean, B, Single_mean, B, Single_to_grouped_mean,B, grouped_to_single_mean,'r', 'LineWidth', 3);
legend(p,'grouped to grouped','single to single','single to grouped','grouped to single')
xlabel('\beta', 'FontSize',30);
ylabel('Transition', 'FontSize',30);
set(gca, 'FontSize',24)
hold off


