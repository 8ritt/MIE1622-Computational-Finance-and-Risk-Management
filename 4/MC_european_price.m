function [callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(S0, K, T, r, mu, sigma, stepNum, numPaths)

% Simulate asset paths for the geometric random walk
S = GRWPaths(S0, mu, sigma, T, stepNum, numPaths);

% Calculate the payoff for each path for a Put
PutPayoffT = max(K-(S(end,:)),0);

% Calculate the payoff for each path for a Call
CallPayoffT = max((S(end,:))-K,0);

% Discount back
putMC_European_Price_1_step = mean(PutPayoffT)*exp(-r*T);
callMC_European_Price_1_step = mean(CallPayoffT)*exp(-r*T);


numSteps = stepNum;

% figure;
% set(gcf, 'color', 'white');
% plot(0:numSteps, S', 'Linewidth', 2);
% title('Geometric Random Walk Paths of the Stock Price', 'FontWeight', 'bold');
% xlabel('Time')
% ylabel('Price')
% 
% 
% figure;
% [frequencyCounts, binLocations] = hist(PutPayoffT, 100);
% bar(binLocations, frequencyCounts,'DisplayName','True Distrubution'); hold on;
% line([mean(PutPayoffT) mean(PutPayoffT)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Mean');
% text(1.2*mean(PutPayoffT), max(frequencyCounts)/1.9, 'Mean Put Payoff');
% title('Put Payoff Distribution')
% xlabel(' Payoff of Put option');
% ylabel('Scenarios');
% 
% figure;
% [frequencyCounts, binLocations] = hist(CallPayoffT, 100);
% bar(binLocations, frequencyCounts,'DisplayName','True Distrubution'); hold on;
% line([mean(CallPayoffT) mean(CallPayoffT)], [0 max(frequencyCounts)/2], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Mean');
% text(1.2*mean(CallPayoffT), max(frequencyCounts)/1.9, 'Mean Call Payoff');
% title('Call Payoff Distribution')
% xlabel(' Payoff of Call option');
% ylabel('Scenarios');
end