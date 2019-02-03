function [daily_value] = daily_portf_value(data_prices, x, days_in_each_period,cash, strategy)
positions_in_years_array =[];
daily_cash_array =[];
for period = 1:length(days_in_each_period)
    
    positions_in_years_array = [positions_in_years_array, repmat(x{strategy,period},1,days_in_each_period(1, period))];
    daily_cash_array =[daily_cash_array,repmat(cash{strategy,period},1,days_in_each_period(1, period))];
end

daily_value = diag(data_prices * positions_in_years_array) + daily_cash_array';

end