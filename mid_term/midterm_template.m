clc;
clear all;
format long

% CSV file with price data
input_file_prices  = 'Daily_closing_prices.csv';

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

% Remove datapoints for year 2014
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Compute means and covariances for Question 2
day_ind_start = 1;
day_ind_end = 39;
cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
mu = mean(cur_returns)';  % Expected returns for Question 2
Q = cov(cur_returns);     % Covariances for Question 2
%%

% Question 1

% Specify quantile level for VaR/CVaR
alf = 0.95;

% Positions in the portfolio
positions = [100 0 0 0 0 0 0 0 200 500 0 0 0 0 0 0 0 0 0 0]';

% Number of assets in universe
Na = size(data_prices,2);

% Number of historical scenarios
Ns = size(data_prices,1);

%%%%% Insert your code here 

stock_loss_1d = - diff(data_prices);

for i = 1:Ns +1 -10
    stock_loss_10d(i,:) = - (data_prices(i+10 -1,:) - data_prices(i,:));
    
end

Ns_10 = Ns + 1 -10;
Ns_1 = Ns - 1;
portf_loss_1d = sort(sum(stock_loss_1d .* positions' , 2));
portf_loss_10d = sort(sum(stock_loss_10d  .* positions', 2));

VaR1  = portf_loss_1d(ceil(Ns_1 * alf));
VaR10 = portf_loss_10d(ceil(Ns_10 * alf));

CVaR1 = (1/(Ns_1*(1-alf))) * ( (ceil(Ns_1*alf)-Ns_1*alf) * VaR1 + sum(portf_loss_1d(ceil(Ns_1*alf)+1:Ns_1)) );
CVaR10 = (1/(Ns_10*(1-alf))) * ( (ceil(Ns_10*alf)-Ns_10*alf) * VaR10 + sum(portf_loss_10d(ceil(Ns_10*alf)+1:Ns_10)) ); 
% 
% VaR1n = mean(portf_loss_1d) + norminv(alf,0,1)*std(portf_loss_1d);
% VaR10n = mean(portf_loss_10d) + norminv(alf,0,1)*std(portf_loss_10d);
% 
% CVaR1n = mean(portf_loss_1d) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_1d);
% CVaR10n = mean(portf_loss_10d) + (normpdf(norminv(alf,0,1))/(1-alf))*std(portf_loss_10d);

day_ind_start_2y = 1;
day_ind_end_2y = 503;
cur_returns_2y = data_prices(day_ind_start_2y+1:day_ind_end_2y,:) ./ data_prices(day_ind_start_2y:day_ind_end_2y-1,:) - 1;
mu_2y = - mean(cur_returns_2y)';  
Q_2y = cov(cur_returns_2y);   

value = data_prices(1,:) * positions;
w_p = (data_prices(1,:) .* positions')' / value;
mu_p = mu_2y' * w_p * value;
std_p = sqrt(w_p' * Q_2y * w_p)  * value;

VaR1n = mu_p + norminv(alf,0,1)*std_p;
CVaR1n = mu_p + (normpdf(norminv(alf,0,1))/(1-alf))*std_p;
%
day_ind_start_10 = 1;
day_ind_end_10 = 495;
cur_returns_10 = data_prices(day_ind_start_2y+9:day_ind_end_2y,:) ./ data_prices(day_ind_start_2y:day_ind_end_2y-9,:) - 1;
mu_10 = - mean(cur_returns_10)';  
Q_10 = cov(cur_returns_10);   

value = data_prices(1,:) * positions;
w_p = (data_prices(1,:) .* positions')' / value;
mu_p_10 = mu_10' * w_p * value;
std_p_10 = sqrt(w_p' * Q_10 * w_p) * value;

VaR10n = mu_p_10 + norminv(alf,0,1)*std_p_10;
CVaR10n = mu_p_10 + (normpdf(norminv(alf,0,1))/(1-alf))*std_p_10;

fprintf('\npart 1\n')
fprintf('Historical 1-day VaR %4.1f%% = $%6.2f,   Historical 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1, 100*alf, CVaR1)
fprintf('    Normal 1-day VaR %4.1f%% = $%6.2f,       Normal 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1n, 100*alf, CVaR1n)
fprintf('Historical 10-day VaR %4.1f%% = $%6.2f,   Historical 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10, 100*alf, CVaR10)
fprintf('    Normal 10-day VaR %4.1f%% = $%6.2f,       Normal 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10n, 100*alf, CVaR10n)


% Plot a histogram of the distribution of losses in portfolio value for 1 day 
figure(1)

[frequencyCounts, binLocations] = hist(portf_loss_1d, 100);
bar(binLocations, frequencyCounts); hold on;
line([VaR1 VaR1], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
%normf = ( 1/(std(portf_loss_1d)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_1d))/std(portf_loss_1d)).^2 ); normf = normf * sum(frequencyCounts)/sum(normf);
%plot(binLocations, normf, 'r', 'LineWidth', 3); hold on;
line([CVaR1 CVaR1], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.');
hold off;
text(1.25*VaR1, max(frequencyCounts)/1.9, 'CVaR1'); text(0.7*VaR1n, max(frequencyCounts)/1.9, 'VaR1');
title('1 day time horizon')
xlabel('Loss')
ylabel('Number of Scenario')

% Plot a histogram of the distribution of losses in portfolio value for 10 days
figure(2)

[frequencyCounts, binLocations] = hist(portf_loss_10d, 100);
bar(binLocations, frequencyCounts); hold on;
line([VaR10 VaR10], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
hold on;
%normf = ( 1/(std(portf_loss_10d)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(portf_loss_10d))/std(portf_loss_10d)).^2 ); normf = normf * sum(frequencyCounts)/sum(normf);
%plot(binLocations, normf, 'r', 'LineWidth', 3); hold on;
line([CVaR10 CVaR10], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-.');
hold off;
text(0.98*CVaR10, max(frequencyCounts)/1.9, 'CVaR10'); text(0.7*VaR10, max(frequencyCounts)/1.9, 'VaR10');
title('10 day time horizon')
xlabel('Loss')
ylabel('Number of Scenario')
%part 2

MSFT_loss = stock_loss_1d(:,1);
AAPL_loss = stock_loss_1d(:,9);
IBM_loss = stock_loss_1d(:,10);

MSFT_loss = sort(MSFT_loss *100);
AAPL_loss = sort(AAPL_loss *200);
IBM_loss = sort(IBM_loss * 500);

VaR1_MSFT  = MSFT_loss(ceil(Ns_1 * alf));
VaR1_AAPL  = AAPL_loss(ceil(Ns_1 * alf));
VaR1_IBM  = IBM_loss(ceil(Ns_1 * alf));
fprintf('\npart 2\n')
fprintf('Historical 1-day MSFT VaR %4.1f%% = $%6.2f\n', 100*alf, VaR1_MSFT)
fprintf('Historical 1-day AAPL VaR %4.1f%% = $%6.2f\n', 100*alf, VaR1_AAPL)
fprintf('Historical 1-day IBM VaR %4.1f%% = $%6.2f\n', 100*alf, VaR1_IBM)


%%
% Question 2

addpath('/Applications/CPLEX_Studio128/cplex/matlab/x86-64_osx');

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025/252;

% Initial portfolio weights
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';
init_value = data_prices(day_ind_end+1,:) * init_positions;
w_init = (data_prices(day_ind_end+1,:) .* init_positions')' / init_value;

% Max Sharpe Ratio portfolio weights
w_Sharpe = [ 0 0 0 0 0 0 0 0.385948690661642 0.172970428625544 0 0 0 0 0 0.003409676869715 0.260942060896445 0 0.185966939781285 0 0]';

% Equal Risk Contribution portfolio weights
w_ERC = [0.049946771209069 0.049951626261681 0.049955739901370 0.049998404150207 0.050000297368719 0.050004255546315 0.050006307026730 0.050007308995726 0.050010525832832 0.050013840015521 0.050014404492514 0.050015932843104 0.050016630302524 0.050017212457105 0.050017600497611 0.050017998351827 0.050018997074443 0.050019598350121 0.050019778113513 0.049946771209069]';

%%%%% Insert your code here 
[info] = efficient_frontier(Na, mu, Q, r_rf, w_init, w_Sharpe, w_ERC);

% Plot for Question 2, Part 1
figure(3)
set(gcf, 'color', 'white');
plot(sqrt(info{1,1}{1,2}), info{1,1}{1,1}, 'k-', 'LineWidth', 1)
hold on;
plot(sqrt(info{1,2}{1,2}), info{1,2}{1,1}, 'rd', 'MarkerSize', 12)
hold on;
plot(sqrt(info{1,3}{1,2}), info{1,3}{1,1}, 'rs', 'MarkerSize', 12)
hold on;
plot(sqrt(info{1,4}{1,2}), info{1,4}{1,1},'mx', 'MarkerSize', 12)
hold on;
plot(sqrt(info{1,5}{1,2}), info{1,5}{1,1}, 'c*','MarkerSize', 12)
hold on;
plot(sqrt(info{1,6}{1,2}), info{1,6}{1,1},'r*', 'MarkerSize', 12)
hold on;
plot(sqrt(info{1,7}{1,2}), info{1,7}{1,1}, 'ro','MarkerSize', 12)
hold on;
plot(sqrt(info{1,8}{1,2}), info{1,8}{1,1},'b^', 'MarkerSize', 12)
hold on;
plot(sqrt(info{1,9}{1,2}), info{1,9}{1,1},'bv', 'MarkerSize', 12)
hold on;
plot(info{1,10}{1,2}, info{1,10}{1,1},'b--', 'LineWidth', 1)
hold off;
%plot(sqrt(diag(Q)), mu, 'b.', 'MarkerSize', 18)
xlabel('2 month Standard deviation');
ylabel('2 month Expected return');
title('Q 2 part 1')
legend('efficient frontier', 'minimum variance portfolio', 'maximum return portfolio', 'equal weight portfolio',...
    'initial portfolio','max Sharpe portfolio','risk free asset', 'equal risk contribution portfolio', ...
    'leveraged equal risk contribution portfolio','efficient frontier with risk free asset')

% Plot for Question 2, Part 2
[info2] = efficient_frontier2(Na, mu, Q);
figure(4);
set(gcf, 'color', 'white');
plot(sqrt(info2{1,1}{1,2}), info2{1,1}{1,1}, 'k-', 'LineWidth', 1)
hold on;
plot(sqrt(info2{1,2}{1,2}), info2{1,2}{1,1}, 'b.', 'MarkerSize', 6)
hold on;
plot(sqrt(diag(Q)), 39 * mu, 'r.', 'MarkerSize', 18)
xlabel('2 month Standard deviation');
ylabel('2 month Expected return');
title('Q 2 Part 2')
legend('efficient frontier', 'random portfolio','individual asset')
hold off;