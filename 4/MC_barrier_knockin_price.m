function [callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = ...
    MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)

% Simulate asset paths for the geometric random walk
S = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);

% S_new = [];
% for j = length(1:numPaths)
%     if sum(S(:,j) >= Sb) ~= 0
%         S_new = [S_new,S(:,j)];
%     end
% end

[S_new] = knock_in(S,Sb);

if isequal(S_new,[])
    PutPayoffT = 0;
    CallPayoffT = 0;
else
    CallPayoffT = max((S_new(end,:))- K,0);
    PutPayoffT = max(K-(S_new(end,:)),0);
end

% Discount back
putMC_Barrier_Knockin_Price_1_step = mean(PutPayoffT)*exp(-r*T);
callMC_Barrier_Knockin_Price_1_step = mean(CallPayoffT)*exp(-r*T);


end