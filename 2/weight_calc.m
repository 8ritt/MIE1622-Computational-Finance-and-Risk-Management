function [w, init_value] = weight_calc(x_init, cur_prices,cash)

init_value = cur_prices * x_init + cash;

% Initial portfolio weights
w = (cur_prices .* x_init')' / init_value;

end