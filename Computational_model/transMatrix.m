function A = transMatrix(domp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the A(i,j) is the number of times that percept i 
% switches to percept j.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pList = domp(:,1); % percept sequences in the time order
numP = length(pList);
A = zeros(4,4); 

for i = 1:numP-1 
    cur_p = pList(i);
    next_p = pList(i+1);
    A(cur_p, next_p) = A(cur_p, next_p)+1;
end

