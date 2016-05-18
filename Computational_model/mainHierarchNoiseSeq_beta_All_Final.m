clear all;
parameters3;
%rand('state', 200);
Dt = 0.1; %ms
tstart = 0;
tend = 1500;
tspan = tstart:Dt:tend;
T = tspan;
N = length(tspan);

nSeq = 500; % number of time series under each given beta value.

b0 = 0.3;
b1 = 0.6;
db = 0.01;
% beta = b0;
params(2) = 0.6; 
params(6) = 1.2;
nb = floor((b1-b0)/db);

%%%%%%%%%%%%%%%%%%%%%
%% Simulating data.%%
%%%%%%%%%%%%%%%%%%%%%
Dom_Single_mean = zeros(nb, 1); % average dominance duration of single eye percepts
Dom_Grouped_mean = zeros(nb, 1); % average dominance duration of single eye percepts
nVisit_Group_mean = zeros(nb, 1); 
Dom_ratiosGroup_mean = zeros(nb, 1); % T_Grouped/ (T_Grouped + T_single)

Dom_Single_std = zeros(nb, 1); % average dominance duration of single eye percepts
Dom_Grouped_std = zeros(nb, 1); % average dominance duration of single eye percepts
nVisit_Group_std = zeros(nb, 1); 
Dom_ratiosGroup_std = zeros(nb, 1); % T_Grouped/ (T_Grouped + T_single)

%%  transition
Single_mean = zeros(nb, 1); % single to single
Grouped_mean = zeros(nb, 1); % grouped to grouped
Single_to_grouped_mean = zeros(nb, 1); % single to grouped
grouped_to_single_mean = zeros(nb, 1); % grouped to single


Single_std = zeros(nb, 1); % single to single
Grouped_std = zeros(nb, 1); % grouped to grouped
Single_to_grouped_std = zeros(nb, 1); % single to grouped
grouped_to_single_std = zeros(nb, 1); % grouped to single


betaSeq = b0:db:b1;

for b=1:nb+1  
    %b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    params(3) = betaSeq(b);
   
    
    dompSingle = zeros(nSeq*50, 2);
    nSingle = 0;
    dompGroup = zeros(nSeq*50, 2);
    nGroup = 0;
    
    nVisRatioGrouped = zeros(nSeq,1); %the ratios of the number of visits to grouped percepts
    TotDomRatioGrouped = zeros(nSeq,1); %ratios of the total dominance of grouped % percepts.
    
    
     trans_ss = zeros(nSeq,1);
    trans_sg = zeros(nSeq,1);
    trans_gg = zeros(nSeq,1);
    trans_gs = zeros(nSeq,1);
    for u = 1:nSeq
        %% Add colored noise;
        %%% Noises for downstrem
        % m  -- number of series noise generated.
        m = 16;
        sig = 0.03; %%% standard deviation.
        noise = ColorNoise(Dt, N, m, sig);
        %%%%%%%%%%%%%%%%%%%
        
        yzero = zeros(16, 1);
        %% %% Euler method
        Y = zeros(N, 16);
        Y(1,:) = yzero';
        for i  = 2: N
            yn = Y(i-1, :)';
            fn = HierarchNoise3(yn, noise(i, :),params);
            ynplus1 = yn + Dt*fn;
            Y(i, :) = ynplus1';
            % T(i) = T(i-1) + Dt;
        end
        %Ys(:, :, u) = Y;
        
        %Analysis
      Y= Y(100:end,:);
         T1 = T(100:end);
        
        
      
        domp = FindDominancePeriods(T1, Y);
   
        
        cut = 3;
        
        domDur = domp(:, 2);
        [r,c] = find(domDur>cut);
        domp = domp(r, :);
        len = length(domp);
        A = transMatrix(domp);
        
        ss = A(1,1)+A(1,2)+A(2,1)+A(2,2);
        sg = A(1,3)+A(1,4)+A(2,3)+A(2,4);
        gs = A(3,1)+A(3,2)+A(4,1)+A(4,2);
        gg = A(3,3)+A(3,4)+A(4,3)+A(4,4);
        
        if (ss+sg >0)
         trans_ss(u) = ss / (ss+ sg);
          trans_sg(u) = 1-  trans_ss(u);
          else 
            trans_ss(u)=0;
            trans_sg(u) = 1;
        end
        if (gg+gs>0)
         trans_gg(u) = gg / (gs + gg);
         trans_gs(u) = 1 - trans_gg(u) ;
        else 
            trans_gg(u)=0;
            trans_gs(u) = 1;
        end
         
         
         perc = domp(:,1);
         [r,c]=find(perc<3);
         singleDom = domp(r,1:2);
         [r,c]=find(perc>2);
         groupedDom = domp(r,1:2);
         
         if ~isempty(groupedDom(:,1))
             numVistGrouped = length(groupedDom(:,1));
             DomGrouped = sum(groupedDom(:,2));
             dompGroup(nGroup+1:nGroup+numVistGrouped,:) = groupedDom;
             nGroup = nGroup + numVistGrouped;
         else
              numVistGrouped = 0;
              DomGrouped = 0;
         end
         
         if ~isempty(singleDom(:,1))
              numVistSingle = length(singleDom(:,1));
              DomSingle = sum(singleDom(:,2));
              dompSingle(nSingle+1:nSingle+numVistSingle,:) = singleDom;
              nSingle = nSingle+numVistSingle;
         else
             numVistSingle = 0;
             DomSingle = 0;
         end
         
        
         nVisRatioGrouped(u) =  numVistGrouped/(numVistGrouped+numVistSingle);
          
         TotDomRatioGrouped(u) = DomGrouped/(DomGrouped+DomSingle);
         
        
    end % End of the nSeq.
    
    dompSingle=dompSingle(1:nSingle,:);
    dompGroup = dompGroup(1:nGroup,:);
    
    % compute mean and variance
    Dom_Single_mean(b) = mean(dompSingle(:,2)); % average dominance duration of single eye percepts
    Dom_Single_std(b) = std(dompSingle(:,2)); % standard deviation of dominance duration of single eye percepts
    
    % compute mean and variance
    %if isempty(dompGroup(:,2)>0)
    Dom_Grouped_mean(b) = mean(dompGroup(:,2)); % average dominance duration of single eye percepts
    Dom_Grouped_std(b) = std(dompGroup(:,2)); % standard deviation of dominance duration of single eye percepts
    
    nVisit_Group_mean(b) = mean(nVisRatioGrouped);
    nVisit_Group_std(b) = std(nVisRatioGrouped);
    
    Dom_ratiosGroup_mean(b) = mean(TotDomRatioGrouped); %
    Dom_ratiosGroup_std(b) = std(TotDomRatioGrouped); %
    Grouped_mean(b) = mean(trans_gg); % grouped to grouped
     grouped_to_single_std(b) = std(trans_gs); % grouped to single
       Grouped_std(b) = std(trans_gg); % grouped to grouped
     grouped_to_single_mean(b) = mean(trans_gs); % grouped to single
%     else
%         
%     end
    Single_mean(b) = mean(trans_ss); % single to single
   
    
    Single_to_grouped_mean(b) = mean(trans_sg); % single to grouped
   
    
    
    Single_std(b) = std(trans_ss); % single to single
  
    Single_to_grouped_std(b) = std(trans_sg); % single to grouped
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end of beta loop


%% Plot the Analysis

B = b0:db:b1;
B = B(1:nb+1);

% Relateion between average dominance duration and coherent coupling
% strength beta.

figure(1)
y1 = Dom_ratiosGroup_mean - Dom_ratiosGroup_std;
y2 = Dom_ratiosGroup_mean + Dom_ratiosGroup_std;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.1 0.1 1]);
hold on
plot(B, Dom_ratiosGroup_mean, 'r', 'LineWidth', 3);
%title('Average Ratios of Grouped percepts', 'FontSize',20)
xlabel('\beta', 'FontSize',30);
ylabel('Predominance', 'FontSize',30);
title('Grouped percepts')
set(gca, 'FontSize',24)
hold off

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





figure(6)
y1 = (Dom_Grouped_mean - Dom_Grouped_std)/10;
y2 = (Dom_Grouped_mean + Dom_Grouped_std)/10;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 1]);
hold on
%plot(B, Dom_Grouped_mean/10, 'r', 'LineWidth', 2,'black')

y1 = (Dom_Single_mean - Dom_Single_std)/10;
y2 = (Dom_Single_mean + Dom_Single_std)/10;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 1]);
hold on


p=plot(B, Dom_Grouped_mean/10, B, Dom_Single_mean/10,'LineWidth', 3);
legend(p,'grouped', 'single eye')
ylabel('Dom. Duration(sec)', 'FontSize',30);
xlabel('\beta', 'FontSize', 30)
set(gca, 'FontSize',24)
hold off






figure(9)

[r,c] = find(~isnan(Grouped_mean));
y1 = Grouped_mean(r) - Grouped_std(r);
y2 = Grouped_mean(r) + Grouped_std(r);
B1 = B(r);
fill_between(B(r), y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.98 0.97 0.97]);
hold on
%plot(B(r), Grouped_mean(r), 'r', 'LineWidth', 3)

[r,c] = find(~isnan(Grouped_mean));
y1 = Single_mean - Single_std;
y2 = Single_mean + Single_std;
y1 = y1(r);
y2 =y2(r);
fill_between(B(r), y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]);
%hold on
y1 = Single_to_grouped_mean - Single_to_grouped_std;
y2 = Single_to_grouped_mean + Single_to_grouped_std;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 1]);
%hold on
plot(B, Single_to_grouped_mean, 'r', 'LineWidth', 3)

%y1 = grouped_to_single_mean - grouped_to_single_std;
%y2 = grouped_to_single_mean + grouped_to_single_std;
fill_between(B, y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.1 0.1 1]);
%hold on
%plot(B, grouped_to_single_mean, 'r', 'LineWidth', 3)
%title('Average Ratios of Grouped percepts', 'FontSize',20)
p=plot(B,Grouped_mean, B, Single_mean, B, Single_to_grouped_mean,B, grouped_to_single_mean,'r', 'LineWidth', 3);
legend(p,'grouped to grouped','single to single','single to grouped','grouped to single')
xlabel('\beta', 'FontSize',30);
ylabel('Transition', 'FontSize',30);
set(gca, 'FontSize',24)
hold off


figure(10)

[r,c] = find(~isnan(Grouped_mean));
y1 = Grouped_mean(r) - Grouped_std(r);
y2 = Grouped_mean(r) + Grouped_std(r);
B1 = B(r);
fill_between(B(r), y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.98 0.97 0.97]);
hold on
%plot(B(r), Grouped_mean(r), 'r', 'LineWidth', 3)

[r,c] = find(~isnan(Grouped_mean));
y1 = Single_mean - Single_std;
y2 = Single_mean + Single_std;
y1 = y1(r);
y2 =y2(r);
fill_between(B(r), y1, y2, 1, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]);

%plot(B, grouped_to_single_mean, 'r', 'LineWidth', 3)
%title('Average Ratios of Grouped percepts', 'FontSize',20)
p=plot(B,Grouped_mean, B, Single_mean, 'r', 'LineWidth', 3);
legend(p,'grouped to grouped','single to single')
xlabel('\beta', 'FontSize',30);
ylabel('Transition', 'FontSize',30);
set(gca, 'FontSize',24)
hold off


%save('result.mat')

