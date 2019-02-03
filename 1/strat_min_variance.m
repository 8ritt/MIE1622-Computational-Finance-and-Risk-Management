function  [x_optimal, cash_optimal] = strat_min_variance(x_init, cash_init, mu, Q, cur_prices)

   n = length(x_init);
   % Optimization problem data
   lb = zeros(n,1);
   ub = inf*ones(n,1);
   A  = ones(1,n);
   b  = 1;
   
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
   
   totalvalue = cur_prices * x_init + cash_init;
   
   value_vector = w_minVar .* totalvalue;
   
   x_optimal = round(value_vector ./ (cur_prices'));
   
   portfolio_changes = x_init - x_optimal;
   
   cash_optimal = cash_init + cur_prices * portfolio_changes;


end