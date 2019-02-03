function  [x_optimal, cash_optimal, borrow] = strat_equal_risk_contr_validation(x_init, cash_init, mu, Q, cur_prices, rf, borrow_init)

% Initial portfolio weights
[w_init, init_value] = weight_calc(x_init,cur_prices,cash_init);
borrow = borrow_init;

global QQQ A_ineq A_eq

QQQ = Q;

%% Add PATH to IPOPT
addpath('/Applications');

% Equality constraints
[m,n] = size(cur_prices);
A_eq = ones(1,n);
b_eq = 1;

% Inequality constraints
A_ineq = [];
b_ineql = [];
b_inequ = [];
           
% Define initial portfolio ("equally weighted" or "1/n portfolio")
w0 = w_init;

options.lb = zeros(1,n);       % lower bounds on variables
options.lu = ones (1,n);       % upper bounds on variables
options.cl = [b_eq' b_ineql']; % lower bounds on constraints
options.cu = [b_eq' b_inequ']; % upper bounds on constraints

% Set the IPOPT options
options.ipopt.jac_c_constant        = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.tol                   = 1e-10;
options.ipopt.print_level           = 0;

% The callback functions
funcs.objective         = @computeObjERC;
funcs.constraints       = @computeConstraints;
funcs.gradient          = @computeGradERC;
funcs.jacobian          = @computeJacobian;
funcs.jacobianstructure = @computeJacobian;

% !!!! Function "computeGradERC" is just the placeholder
% !!!! You need to compute the gradient yourself
  
%% Run IPOPT
[wsol, info] = ipopt(w0',funcs,options);

% Make solution a column vector
if(size(wsol,1)==1)
    w_erc = wsol';
else
    w_erc = wsol;
end

value_vector = w_erc * init_value;

x_optimal = round(value_vector ./ cur_prices');

x_change = x_init - x_optimal;

balancing_cost = cur_prices * x_change;

cash_optimal = cash_init + balancing_cost;

% % Compute return, variance and risk contribution for the ERC portfolio
ret_ERC = dot(mu, w_erc);
var_ERC = w_erc'*Q*w_erc;
RC_ERC = (w_erc .* ( Q*w_erc )) / sqrt(w_erc'*Q*w_erc);
% 



RC_init  = (w0 .* ( Q*w0 )) / sqrt(w0' * Q * w0);
% 


fprintf('\n\nAsset risk contributions for ERC:\n')
RC_ERC
% 

