clc;
clear all;
format long

addpath('/Applications');
addpath('/Applications/CPLEX_Studio128/cplex/matlab/x86-64_osx');
% Input files
input_file_prices  = 'Daily_closing_prices20082009.csv';

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     % Read time periods
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
dates_array = dates_array(:,1:3);
%%
% Find the number of trading days in Nov-Dec 2014 and
% compute expected return and covariance matrix for period 1
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2007));
cur_returns0 = data_prices(day_ind_start0+1:day_ind_end0,:) ./ data_prices(day_ind_start0:day_ind_end0-1,:) - 1;
mu = mean(cur_returns0)';
Q = cov(cur_returns0);

% Remove datapoints for year 2014
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Initial positions in the portfolio
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';
%%
% Initial value of the portfolio
init_value = data_prices(1,:) * init_positions;
fprintf('\nInitial portfolio value = $ %10.2f\n\n', init_value);

% Initial portfolio weights
w_init = (data_prices(1,:) .* init_positions')' / init_value;

% Number of periods, assets, trading days
N_periods = 6*length(unique(dates_array(:,1))); % 6 periods per year
N = length(tickers);
N_days = length(dates);

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;
% Annual risk-free rate for years 2008-2009 is 4.5%
r_rf2008_2009 = 0.045;

rf = r_rf2008_2009;
% initial leverage used

borrow_last_period = 0;

% Number of strategies
strategy_functions = {'strat_buy_and_hold' 'strat_equally_weighted' 'strat_min_variance' 'strat_max_Sharpe' 'strat_equal_risk_contr' 'strat_lever_equal_risk_contr' 'strat_robust_optim'};
strategy_names     = {'Buy and Hold' 'Equally Weighted Portfolio' 'Minimum Variance Portfolio' 'Maximum Sharpe Ratio Portfolio' 'Equal Risk Contributions Portfolio' 'Leveraged Equal Risk Contributions Portfolio' 'Robust Optimization Portfolio'};
%N_strat = 5; % comment this in your code
N_strat = length(strategy_functions); % uncomment this in your code
fh_array = cellfun(@str2func, strategy_functions, 'UniformOutput', false);


for (period = 1:N_periods)
   % Compute current year and month, first and last day of the period
   if(dates_array(1,1)==08)
       cur_year  = 8 + floor(period/7);
   else
       cur_year  = 2008 + floor(period/7);
   end
   cur_month = 2*rem(period-1,6) + 1;
   day_ind_start = find(dates_array(:,1)==cur_year & dates_array(:,2)==cur_month, 1, 'first');
   day_ind_end = find(dates_array(:,1)==cur_year & dates_array(:,2)==(cur_month+1), 1, 'last');
   fprintf('\nPeriod %d: start date %s, end date %s\n', period, char(dates(day_ind_start)), char(dates(day_ind_end)));

   % Prices for the current day
   cur_prices = data_prices(day_ind_start,:);

   % Execute portfolio selection strategies
   for(strategy = 1:N_strat)

      % Get current portfolio positions
      if(period==1)
         curr_positions = init_positions;
         curr_cash = 0;
         portf_value{strategy} = zeros(N_days,1);
      else
         curr_positions = x{strategy,period-1};
         curr_cash = cash{strategy,period-1};
      end
      
      if period == 1
          borrow_last_period = 0;
      else
          borrow_last_period = borrow{strategy,period-1};
      end
      % Compute strategy
      [x{strategy,period} cash{strategy,period}, borrow{strategy,period}] = fh_array{strategy}(curr_positions, curr_cash, mu, Q, cur_prices, rf, borrow_last_period);
      
      borrow_last_period = borrow;
      % Verify that strategy is feasible (you have enough budget to re-balance portfolio)
      % Check that cash account is >= 0
      % Check that we can buy new portfolio subject to transaction costs

      %%%%%%%%%%% Insert your code here %%%%%%%%%%%%
      rebalancing_transaction = curr_positions - x{strategy, period};
      transaction_cost = cur_prices * abs(rebalancing_transaction) * 0.005;
      cash{strategy, period} = cash{strategy, period} - transaction_cost;
      
      if cash{strategy, period} < 0
          [M, I] = max(x{strategy, period});
          stock_to_sell_index = I;
          stock_to_sell_price = cur_prices(1, stock_to_sell_index);
          number_of_shares_sell = ceil(abs(cash{strategy, period})/ (1 - 0.005) / stock_to_sell_price);
          
          new_cash = number_of_shares_sell * stock_to_sell_price + cash{strategy, period};
          cash{strategy, period} = round(new_cash, 2);
          x_new = x{strategy, period}(stock_to_sell_index , 1) - number_of_shares_sell;
          
          x{strategy, period}(stock_to_sell_index , 1) = x_new;
      end
      % Compute portfolio value
      portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period} - borrow{strategy,period} *(1+rf/6);

      fprintf('   Strategy "%s", value begin = $ %10.2f, value end = $ %10.2f\n', char(strategy_names{strategy}), portf_value{strategy}(day_ind_start), portf_value{strategy}(day_ind_end));

   end
      
   % Compute expected returns and covariances for the next period
   cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
   mu = mean(cur_returns)';
   Q = cov(cur_returns);

end
%%
% Plot results
% figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%
figure(1);
days_in_each_period = [];
for year = [2008, 2009]
    for month = 1:12
        days_in_month = length(find(dates_array(:,1)==year & dates_array(:,2)==month));
        days_in_each_period = [days_in_each_period days_in_month];
    end
end

days_in_each_period = sum(reshape(days_in_each_period, 2, 12));

days =1:length(dates);

hold on
plot(days, portf_value{1},'DisplayName', 'strat buy and hold')
plot(days, portf_value{2},'DisplayName', 'strat equally weighted')
plot(days, portf_value{3}, 'DisplayName', 'strat min variance')
plot(days, portf_value{4}, 'DisplayName', 'strat max Sharpe')
plot(days, portf_value{5}, 'DisplayName', 'strat equal risk contr')
plot(days, portf_value{6}, 'DisplayName', 'strat leverage equal risk contri')
plot(days, portf_value{7}, 'DisplayName', 'strat robust optim')

legend('show')
title('Daily Value of Portfolio of Each Strategy');
xlabel('Number of Days');
ylabel('Value of Portfolio');
hold off



figure(2)
[weight_array_strat3] = weight_calculation(x,data_prices,days_in_each_period,3,cash);
weight_array_strat3 = [w_init,weight_array_strat3];
plot(0:12,weight_array_strat3)
legend(headers{1,1}(2:21,1))
title('Portfolio Allocations of strat min vairance')
xlabel('Period')
ylabel('Weight of asset')


figure(3)
[weight_array_strat4] = weight_calculation(x,data_prices,days_in_each_period,4,cash);
weight_array_strat4 = [w_init,weight_array_strat4];
plot(0:12,weight_array_strat4)
legend(headers{1,1}(2:21,1))
title('Portfolio Allocations of strat max Sharpe')
xlabel('Period')
ylabel('Weight of asset')

figure(4)
[weight_array_strat7] = weight_calculation(x,data_prices,days_in_each_period,7,cash);
weight_array_strat7 = [w_init,weight_array_strat7];
plot(0:12,weight_array_strat7)
legend(headers{1,1}(2:21,1))
title('Portfolio Allocations of strat robust optim')
xlabel('Period')
ylabel('Weight of asset')
