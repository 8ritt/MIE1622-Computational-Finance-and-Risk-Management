function  [x_optimal, cash_optimal, borrow] = strat_lever_equal_risk_contr(x_init, cash_init, mu, Q, cur_prices,rf, borrow_init)

payment = borrow_init * (1+rf/6)*(1+0.005); % payback to loan plus interest plus transaction fee

% sell stocks to get payback money
[w_init, init_value] = weight_calc(x_init,cur_prices,cash_init);

value_sell_vector = w_init * payment;
shares_sell_vector = ceil(value_sell_vector ./ (cur_prices'));
money_sold = cur_prices * shares_sell_vector * (1-0.005);

% borrow equal amount of money of port value after payback
cash_post_payback = cash_init + money_sold - borrow_init * (1 + rf/6);
x_post_payback = x_init - shares_sell_vector;
[w_init, port_value] = weight_calc(x_post_payback,cur_prices,cash_post_payback);

borrow = port_value;
cash_pre_balancing = borrow + cash_post_payback;

[x_post_balancing, cash_post_balancing, borrow_optimal] = strat_equal_risk_contr(x_post_payback, cash_pre_balancing, mu, Q, cur_prices, rf, borrow);

x_optimal = x_post_balancing;
cash_optimal = cash_post_balancing;

end


