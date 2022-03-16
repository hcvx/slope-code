function output = generate_lambda(n, opt, para)
% generate_lambda: generate regularizing sequence

switch opt
    case 0
        % lasso
        lambda_discrete = ones(1,n);
    case 1
        % discrete
        lambda_finite = para.lambda_finite;
        lambda_finite_prob = para.lambda_finite_prob;
        lambda_finite_cumsum = [0 cumsum(lambda_finite_prob)];
        num_of_lambda_finite = length(lambda_finite);
        lambda_discrete = ones(1,n);
        for i=1:num_of_lambda_finite
            lambda_discrete(fix(lambda_finite_cumsum(i)*n)+1:fix(lambda_finite_cumsum(i+1)*n))=lambda_finite(i);
        end
    case 2
        % uniform
        lambda_discrete = 1*(1:1:n)/n;
    case 3
        % BHq
        q = para.q;
        lambda_discrete = 1*icdf('Normal',1-q*(1:1:n)/(2*n),0,1);
    case 4
        % optimal design
        p_array = linspace(0,1,n);
        absy_array = para.absy_array;
        u_array = para.u_array;
        tau1 = para.tau1;
        lambda_discrete = optlambda_func(absy_array, u_array, p_array, para)/tau1;
end
[lambda_discrete,~] = sort(lambda_discrete,'descend');
output = lambda_discrete;
end