function [info] = efficient_frontier2(Na, mu, Q)
mu = mu * 39;
n = Na;
% Optimization problem data
lb = zeros(n,1);
ub = inf*ones(n,1);
A  = ones(1,n);
b  = 1;
%% Minimum Variance 
% Compute minimum variance portfolio
cplex1 = Cplex('min_Variance');
cplex1.addCols(zeros(n,1), [], lb, ub);
cplex1.addRows(b, A, b);
cplex1.Model.Q = 2*Q;
cplex1.Param.qpmethod.Cur = 6; % concurrent algorithm
cplex1.Param.barrier.crossover.Cur = 1; % enable crossover
cplex1.DisplayFunc = []; % disable output to screen
cplex1.solve();

% Display minimum variance portfolio
w_minVar = cplex1.Solution.x;
var_minVar = w_minVar' * Q * w_minVar;
ret_minVar = mu' * w_minVar;

%% Maximum return 
% Compute maximum return portfolio
cplex2 = Cplex('max_Return');
cplex2.Model.sense = 'maximize';
cplex2.addCols(mu, [], lb, ub);
cplex2.addRows(b, A, b);
cplex2.Param.lpmethod.Cur = 6; % concurrent algorithm
cplex2.Param.barrier.crossover.Cur = 1; % enable crossover
cplex2.DisplayFunc = []; % disable output to screen
cplex2.solve();

% Display maximum return portfolio
w_maxRet = cplex2.Solution.x;
var_maxRet = w_maxRet' * Q * w_maxRet;
ret_maxRet = mu' * w_maxRet;

% Target returns
targetRet = linspace(ret_minVar,ret_maxRet,20);
%% Efficient Frontier no short 
% Compute efficient frontier
cplex3 = Cplex('Efficient_Frontier');
cplex3.addCols(zeros(n,1), [], lb, ub);
cplex3.addRows(targetRet(1), mu', inf);
cplex3.addRows(b, A, b);
cplex3.Model.Q = 2*Q;
cplex3.Param.qpmethod.Cur = 6; % concurrent algorithm
cplex3.Param.barrier.crossover.Cur = 1; % enable crossover
cplex3.DisplayFunc = []; % disable output to screen

w_front = [];
for i=1:length(targetRet)
    cplex3.Model.lhs(1) = targetRet(i);
    cplex3.solve();
    w_front = [w_front cplex3.Solution.x];
    var_front(i) = w_front(:,i)' * Q * w_front(:,i);
    ret_front(i) = mu' * w_front(:,i);
end

%% random portfolios
rnd_portf = zeros(n,1000);

for j = 1:1000
    rnd_portf(:,j) = rand(20,1);
end

w_rnd_portf = rnd_portf ./ (sum(rnd_portf));

for k=1:length(w_rnd_portf)
    w_rnd = w_rnd_portf(:,k);
    var_rnd(k) = w_rnd' * Q * w_rnd;
    ret_rnd(k) = mu' * w_rnd;
end



%%
info ={};
info{1,1} = {ret_front, var_front};
info{1,2} = {ret_rnd, var_rnd};
end

