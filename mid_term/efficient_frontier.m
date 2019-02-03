function [info] = efficient_frontier(Na, mu, Q, rf, w_init, w_Sharpe, w_ERC)

% data for 20 stocks
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

%% equally weighted
w_EQW = 1/n*(ones(n,1));
% 1/n portfolio return
ret_equal_weight = dot(mu, w_EQW);
% 1/n portfolio variance
var_equal_weight = w_EQW' * Q * w_EQW;


%% intial portfolio
ret_init_portf = dot(mu, w_init);
var_init_portf = w_init' * Q * w_init;


%% Maximum Sharpe
ret_Sharpe = dot(mu, w_Sharpe);
var_Sharpe = w_Sharpe' * Q * w_Sharpe;


%% Risk Free
ret_rf = rf ;
var_rf = 0;


%% Equal Risk Contributor
ret_ERC = dot(mu, w_ERC);
var_ERC = w_ERC' * Q * w_ERC;


%% Leverage ERC

ret_LERC = 2 * ret_ERC - ret_rf;
var_LERC = 2^2 * var_ERC;


%% efficient frontier shorting allowed
ret_frontRf(1) = ret_rf;
var_frontRf(1) = var_rf;

ret_frontRf(2) = ret_Sharpe;
var_frontRf(2) = sqrt(var_Sharpe);

ret_frontRf(3) = 2 *(ret_Sharpe - ret_rf);
var_frontRf(3) = 2 * (sqrt(var_Sharpe)-sqrt(var_rf));
% for j = 2:20
%     j = j-1;
%     CML_slope = (ret_Sharpe - rf)/(var_Sharpe - var_rf);
%     ret_frontRf(j) = CML_slope * (j/10) * (var_Sharpe - var_rf);
%     var_frontRf(j) = (j/10)^2 * (var_Sharpe - var_rf);
% end

%%
info ={};
info{1,1} = {ret_front, var_front};
info{1,2} = {ret_minVar, var_minVar};
info{1,3} = {ret_maxRet, var_maxRet};
info{1,4} = {ret_equal_weight, var_equal_weight};
info{1,5} = {ret_init_portf, var_init_portf};
info{1,6} = {ret_Sharpe, var_Sharpe};
info{1,7} = {ret_rf, var_rf};
info{1,8} = {ret_ERC, var_ERC};
info{1,9} = {ret_LERC, var_LERC};
info{1,10} = {ret_frontRf, var_frontRf};



end

