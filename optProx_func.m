function [absy_array, u_array, f_eval, Expe_deri] = optProx_func(num_of_grid, para)

sparsity = para.sparsity;
delta = para.delta;

cond_mean_func = @(y) condX_mean_func(y, para);
pdf_y_func = @(y) pdf_of_X(y, para);
pdf_absy_func = @(y) pdf_y_func(y) + pdf_y_func(-y);

% find grid point

max_y = invCDF_of_absX(1-10^-6, pdf_absy_func);

abs_yk = linspace(0, max_y, num_of_grid+1);


abs_yk = abs_yk(1:num_of_grid);
Delta_k = abs_yk - [0 abs_yk(1:num_of_grid-1)];
abs_yk_center = (abs_yk+[0 abs_yk(1:num_of_grid-1)])/2;
% 2. calculate PDF of |y_k|

pdf_yk_plus = pdf_y_func(abs_yk_center);
pdf_yk_minus = pdf_y_func(-abs_yk_center);
pdf_abs_yk = pdf_yk_plus+pdf_yk_minus;

% calculate the probability of each bin
prob_yk_plus = pdf_yk_plus.*Delta_k;
prob_yk_minus = pdf_yk_minus.*Delta_k;

% 3. calculate conditional expectation of y_k and -y_k
cond_mean_plus = cond_mean_func(abs_yk);
cond_mean_minus = cond_mean_func(-abs_yk);

% 4. calculate parameters for QP
A1 = diag(ones(1,num_of_grid))-diag(ones(1,num_of_grid-1),-1);
A2 = -diag(ones(1,num_of_grid))+diag(ones(1,num_of_grid-1),-1);
A3 = pdf_abs_yk-[pdf_abs_yk(2:num_of_grid) 0];
b1 = Delta_k';
b2 = zeros(num_of_grid,1);
b3 = delta;
h = (prob_yk_plus+prob_yk_minus)';
f = ((cond_mean_minus.*prob_yk_minus - cond_mean_plus.*prob_yk_plus))';
u0 = sparsity*abs_yk'.*rand(size(abs_yk'));
lb = zeros(num_of_grid,1);
ub = max(abs_yk)*ones(num_of_grid,1);

%% optimization
Aeq = [];
beq = [];
A = [A1;A2;A3];
b = [b1;b2;b3];

options = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex', 'MaxIter', 1*10^3, 'TolFun', 10^-15, 'TolX', 10^-15);

[u_opt, f_eval] = quadprog(diag(h),f',A,b,Aeq,beq,lb,ub,u0,options);

f_eval = 2*f_eval;

absy_array = abs_yk;
u_array = u_opt';
Expe_deri = A3*u_opt;
end