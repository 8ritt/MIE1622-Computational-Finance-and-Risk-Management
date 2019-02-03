function [weight_array] = weight_calculation(x,data_prices,days_in_each_period,strategy,cash)

weight_array =zeros(20,12);
j = 1;
for i = 1:12
    shares_array = x{strategy,i};
    cur_prices = data_prices(j,:);
    j = j + days_in_each_period(1,i);
    total_value = cur_prices * shares_array + cash{strategy,i};
    value_vector = diag(shares_array * cur_prices);
    weight_array(:,i) = value_vector ./ total_value;
end