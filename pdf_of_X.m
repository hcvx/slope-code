function fx = pdf_of_X(x, para)

switch para.X_sgldist_name
    case 'gauss'
        sparsity = para.sparsity;
        sigma_H = para.sigma_H;
        sigma_X0 = para.sigma_X0;

        sigma_H0 = sigma_H;
        sigma_H1 = (sigma_H^2+sigma_X0^2)^0.5;
        mu_H0 = 0;
        mu_H1 = para.mu_X0;
        fx = sparsity * normpdf(x,mu_H1,sigma_H1) + (1-sparsity) * normpdf(x,mu_H0,sigma_H0);
        
    case 'discrete'
        sparsity = para.sparsity;
        sigma_H = para.sigma_H;
        sigma_H0 = sigma_H;
        sigma_H1 = sigma_H;
        mu_H0 = 0;
        beta_discrete = para.beta_discrete;
        prob_discrete = para.prob_discrete;
        k = numel(prob_discrete);
        fx = 0;
        for i=1:k
            fx = fx + sparsity * prob_discrete(i) * normpdf(x, beta_discrete(i), sigma_H1);
        end
        fx = fx + (1-sparsity) * normpdf(x, mu_H0, sigma_H0);
        
    case 'uniform'
        sparsity = para.sparsity;
        sigma_H = para.sigma_H;
        a = para.a;
        b = para.b;
        
        fx_zero = normpdf(x, 0, sigma_H);
        fx_not_zero = (normcdf((b - x)/sigma_H) - normcdf((a - x)/sigma_H)) / (b - a);
        fx = (1-sparsity) * fx_zero + sparsity * fx_not_zero;
end


end