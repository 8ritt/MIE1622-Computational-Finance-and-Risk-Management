function [call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma)

d1 = 1/(sigma * sqrt(T)) * (log(S0/K) + (r + (sigma^2) /2) * (T));
d2 = d1 - sigma * sqrt(T);

call_BS_European_Price = normcdf(d1)*S0 - normcdf(d2) * K * exp(-r*(T));
putBS_European_Price = normcdf(-d2) * K * exp(-r *(T)) - normcdf(-d1)*S0;

end
