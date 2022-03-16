function lambda_array = optlambda_func(absy_array, u_array, p_array, para)

% optlambda_func is used to sample lambda array based on the input optimal
% design [absy_array, u_array] obtained from optProx_func
% p_array:  input CDF
% lambda_array: output quantile

list_dim = 10^4; 

para1 = para;
para1.sigma_H = para.design_sigma_z;

input_absy_array = inverse_CDF_absy(p_array, list_dim, para1);
output_u_array = opt_u_func(absy_array, u_array, input_absy_array);

num_of_u_zero = sum(output_u_array==0);
lambda_array_ori = input_absy_array - output_u_array;
min_lambda = lambda_array_ori(min(length(lambda_array_ori), num_of_u_zero+1));

lambda_array = input_absy_array - output_u_array;
lambda_array = max(lambda_array, min_lambda);
end
