function  [x_optimal, cash_optimal, borrow] = strat_equally_weighted(x_init, cash_init, mu, Q, cur_prices, rf, borrow_init)  
    
    n = length(x_init);
    weight = 1/n;
    totalvalue = cur_prices * x_init + cash_init;
    
    valuei = weight * totalvalue;
    value_vector = valuei * ones(n,1);
    %value_vector(:,:) = valuei;
    
    x_optimal = round(value_vector ./ (cur_prices'));
    
    portfolio_changes = x_init - x_optimal;
    
    w_ew = weight * ones(n, 1);
    cash_optimal = cash_init + cur_prices * portfolio_changes;
    
    borrow = borrow_init;
    
end