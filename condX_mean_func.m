function condx = condX_mean_func(x, para)
% compute E(beta|X), where X = beta + sigma_H*w, w~N(0,1), using Tweedie's
% formula:
%   E(beta|X) = X + sigma_H * P'(X) / P(X)
switch para.X_sgldist_name
    case 'gauss'
        sparsity = para.sparsity;
        mu_S = para.mu_X0;
        sigma_X0 = para.sigma_X0;
        sigma_H = para.sigma_H;
        condx = sparsity*(mu_S+sigma_X0^2/(sigma_X0^2+sigma_H^2)*(x-mu_S))./(sparsity + (1-sparsity)*(sigma_H^2+sigma_X0^2)^0.5/sigma_H*exp(0.5*(x-mu_S).^2/(sigma_H^2+sigma_X0^2)-0.5*x.^2/sigma_H^2));
    case 'discrete' 
        sparsity = para.sparsity;
        sigma_H = para.sigma_H;

        beta_discrete = para.beta_discrete;
        prob_discrete = para.prob_discrete;
        k = numel(prob_discrete);
        
        
%         deri_pdfx = 0;
%         for i=1:k
%             deri_pdfx = deri_pdfx + sparsity * prob_discrete(i) * normpdf(x, beta_discrete(i), sigma_H) .* (- (x - beta_discrete(i)) / sigma_H^2 );
%         end
%         deri_pdfx = deri_pdfx + (1 - sparsity) * normpdf(x, 0, sigma_H) .* (- x / sigma_H^2 );
%         pdfx = pdf_of_X(x, para);
%         condx = x + sigma_H^2 * deri_pdfx ./ pdfx;
        
        sum_1 = zeros(size(x));
        for i=1:k
            sum_1 = sum_1 + prob_discrete(i) * exp(- (beta_discrete(i)^2 - 2*beta_discrete(i)*x) / (2*sigma_H^2) );
        end        
        
        condx = zeros(size(x));
        for i=1:k
            sum_2 = zeros(size(x));
            for j=1:k
                if j~=i
                    sum_2 = sum_2 + prob_discrete(j) * exp( - ( (beta_discrete(j)^2 - beta_discrete(i)^2)  - 2*(beta_discrete(j) - beta_discrete(i)) * x) / (2*sigma_H^2) );
                else
                    sum_2 = sum_2 + prob_discrete(i);
                end
            end
            if sparsity == 1
                condx = condx + prob_discrete(i) * beta_discrete(i) * (sparsity*sum_2).^-1;
            else
                condx = condx + prob_discrete(i) * beta_discrete(i) * ( (1-sparsity) * exp(( (beta_discrete(i)^2)  - 2*(beta_discrete(i)) * x) / (2*sigma_H^2) ) + sparsity*sum_2).^-1;
            end
            
            %condx1 = condx1 + prob_discrete(i) * beta_discrete(i) ./ ( exp(( (beta_discrete(i)^2)  - 2*(beta_discrete(i)) * x) / (2*sigma_H^2) ) .* ( (1-sparsity) + sparsity * sum_1  ));
        end
        condx = condx * sparsity;
    case 'uniform'
        sigma_H = para.sigma_H;
        sparsity = para.sparsity;
        a = para.a;
        b = para.b;
        
        deri_pdfx = (1-sparsity)*  normpdf(x, 0, sigma_H) .* (- x / sigma_H^2 )  + sparsity * 1/((b-a)*sigma_H) * ( normpdf((a - x)/sigma_H) - normpdf((b - x)/sigma_H) );
        pdfx = pdf_of_X(x, para);
        condx = x + sigma_H^2 * deri_pdfx ./ pdfx;
end

end