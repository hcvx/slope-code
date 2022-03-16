clear;

SNR = 0.1;

n = 1000;

delta = 0.5;

m = fix(delta*n);

sparsity = 0.2;
para.sparsity=sparsity;
k_sparsity = fix(n*sparsity);

% choose signal prior
X_sgldist_name = 'gauss';
switch X_sgldist_name
    case 'gauss'
        para.mu_X0 = 1;
        para.sigma_X0 = 0.0;
        para.sgl_sq = para.sparsity*(para.sigma_X0^2+para.mu_X0^2);
    case 'discrete'       
        para.beta_discrete = linspace(0, 5, 10);
        para.prob_discrete = 1/numel(para.beta_discrete) * ones(1, numel(para.beta_discrete));
        para.sgl_sq = para.sparsity * para.beta_discrete.^2 * para.prob_discrete';
end
para.X_sgldist_name = X_sgldist_name;

sigma = ( para.sgl_sq / SNR )^0.5;

q = 0.55;
num_of_lambda = n;
lambda_opt = 4;

para.q = q;
para.n = n;
para.delta=delta;
para.lambda_opt = lambda_opt;

lambda_finite = [1 2];
lambda_finite_prob = [0.5 0.5];
lambda_finite_cumsum = [0 cumsum(lambda_finite_prob)];

para.lambda_finite = lambda_finite;
para.lambda_finite_prob = lambda_finite_prob;
para.CDF_lambda_finite = [lambda_finite_cumsum 1];

para.sigma=sigma;

%% generate optimal regularizing sequence
num_of_grid = 2048;
para.num_of_grid =num_of_grid;
lambda_strength = 1;

% find the optimal one-dimensional denoiser
[absy_array, u_array, sigma_z, tau1] = optDesign_func(num_of_grid, para);

para.absy_array = absy_array;
para.u_array = u_array;
para.tau1 = tau1;
para.design_sigma_z = sigma_z;

% generate optimal regularizing sequence
lambda_discrete = generate_lambda(num_of_lambda, para.lambda_opt, para);
lambda_opt = lambda_strength*lambda_discrete;

% calculate theoretical minimum MSE
MSE_design_theo = (sigma_z^2-sigma^2)*delta;
%% Monte carlo simulation
Iter = 10;
MSE_Iter = zeros(1,Iter);

for iter = 1:Iter

z = sigma*randn(m,1);
z1 = randn(m,1);

randomMask = zeros(n,1);
randomMask(randperm(n,k_sparsity))=1;


switch para.X_sgldist_name
    case 'gauss'
        x0 = (para.sigma_X0*randn(n,1)+para.mu_X0*ones(n,1)).*randomMask;
    case 'discrete'
        prob_discrete_CDF = cumsum(para.prob_discrete);
        ind_mtx = rand(n, 1) * ones(1, numel(para.prob_discrete)) > ones(n, 1) * prob_discrete_CDF;
        idx = sum(ind_mtx, 2) + 1;
        x0 = para.beta_discrete(idx)';
        x0 = x0 .* randomMask;
end

A = randn(m,n)/m^0.5;
y = A*x0+z;

[x1,~] = Adlas(A,y,lambda_opt);
MSE_Iter(iter) = norm(x1-x0)^2/n;
end
mean_MSE_Iter = mean(MSE_Iter);

diff = MSE_design_theo - mean_MSE_Iter;



