function  [x_optimal, cash_optimal, borrow] = strat_max_Sharpe(x_init, cash_init, mu, Q, cur_prices, rf, borrow_init)
    
    borrow = borrow_init;
    n = length(x_init);
    nn = n + 1;
    
    QQ = zeros(nn,nn);
    QQ(1:n,1:n) = 2 * Q;
    mu = mu';
    
    rf = rf / 252;
    A = zeros(2 , nn);
    A(1, 1:n) = mu - rf;
    
    % if all mu-rf elements are negative
    if sum((mu-rf) <= 0) == n
        % hold portfolio position dont rebalance
        x_optimal = x_init;
        cash_optimal = cash_init;
        
    else
        
        A(2, 1:n) =  1;
        A(2, nn) = - 1;
        b = [1; 0];
        lb = zeros(nn, 1);
        ub = inf * ones(nn, 1);
        
        
        cplex1 = Cplex('max_sharpe_ratio');
        cplex1.addCols(zeros(nn,1), [], lb, ub);
        cplex1.addRows(b, A, b);
        cplex1.Model.Q = QQ;
        cplex1.Param.qpmethod.Cur = 6; % concurrent algorithm
        cplex1.Param.barrier.crossover.Cur = 1; % enable crossover
        cplex1.DisplayFunc = []; % disable output to screen
        cplex1.solve();
        
        solution = cplex1.Solution.x;
        
        w_maxSharpe = solution(1:n, 1);
        
        w_maxSharpe = w_maxSharpe ./ solution(nn, 1);
        
        totalvalue = cur_prices * x_init + cash_init;
        value_vector = w_maxSharpe .* totalvalue;
        
        x_optimal = round(value_vector ./ (cur_prices'));
        
        portfolio_changes = x_init - x_optimal;
        
        cash_optimal = cash_init + cur_prices * portfolio_changes;
        
    end
end