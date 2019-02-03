function [S_new] = knock_in(S,threshold)
[m,n] = size(S);
S_new = [];

for i = 1:n
    if sum(S(:,i) >= threshold) > 0
        S_new(:,end+1) = S(:,i);
    end
end
end