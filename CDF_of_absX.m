function Fx = CDF_of_absX(x, para)


switch para.X_sgldist_name
    case 'gauss'
        sparsity = para.sparsity;
        sigma_H = para.sigma_H;

        sigma_H0 = sigma_H;
        sigma_H1 = (sigma_H^2+para.sigma_X0^2)^0.5;
        mu_H0 = 0;
        mu_H1 = para.mu_X0;
        Fx = sparsity*(normcdf(abs(x),mu_H1,sigma_H1)-normcdf(-abs(x),mu_H1,sigma_H1))+(1-sparsity)*(normcdf(abs(x),mu_H0,sigma_H0)-normcdf(-abs(x)-10^-8,mu_H0,sigma_H0));
    case 'discrete'
        sparsity = para.sparsity;
        sigma_H = para.sigma_H;
        sigma_H0 = sigma_H;
        sigma_H1 = sigma_H;
        mu_H0 = 0;
        beta_discrete = para.beta_discrete;
        prob_discrete = para.prob_discrete;
        k = numel(prob_discrete);
        Fx = 0;
        for i=1:k
            Fx = Fx + sparsity * prob_discrete(i) * ( normcdf(abs(x), beta_discrete(i), sigma_H1) - normcdf(-abs(x), beta_discrete(i), sigma_H1) );
        end
        Fx = Fx + (1-sparsity) * (normcdf(abs(x),mu_H0, sigma_H0)-normcdf(-abs(x),mu_H0, sigma_H0));
    case 'uniform'   
        pdf_y_func = @(y) pdf_of_X(y, para);
        pdf_absy_func = @(y) pdf_y_func(y) + pdf_y_func(-y);  
        Fx = zeros(size(x));
        for j=1:numel(x);
            Fx(j) = integral(pdf_absy_func, 0, abs(x(j)), 'RelTol', 1e-10, 'Waypoints', 0);
        end
end


end