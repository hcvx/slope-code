function x_ary = invCDF_of_X(CDF_ary, pdf_func)
% find the quantile for X, the support of X is R
% pdf_X_func: pdf of the distribution
% CDF_ary: input CDFs
% x_ary: output quantiles

x_ary = zeros(size(CDF_ary));
for i=1:numel(CDF_ary)
    if CDF_ary(i) == 1
        x_ary(i) = inf;
        continue;
    end
    
    if CDF_ary(i) == 0
        x_ary(i) = 0;
        continue;
    end
    
    % compute P(X<=0)
    CDF_0 = integral(@(y) pdf_func(y), -inf, 0, 'RelTol', 1e-20);
    
    if CDF_ary(i) < CDF_0
        min_X = -2;
        max_X = 0;
        CDF_min_X = integral(@(y) pdf_func(y), -inf, min_X, 'RelTol', 1e-20);
        while CDF_ary(i) < CDF_min_X
            max_X = min_X;
            min_X = min_X * 1.5;
            CDF_min_X = integral(@(y) pdf_func(y), -inf, min_X, 'RelTol', 1e-20);
        end        
    else
        min_X = 0;
        max_X = 2;
        CDF_max_X = integral(@(y) pdf_func(y), -inf, max_X, 'RelTol', 1e-20);
        while CDF_ary(i) > CDF_max_X
            min_X = max_X;
            max_X = max_X * 1.5;
            CDF_max_X = integral(@(y) pdf_func(y), -inf, max_X, 'RelTol', 1e-20);
        end        
    end   
    
    cur_X = (min_X + max_X) / 2;
    cur_CDF = integral(@(y) pdf_func(y), -inf, cur_X, 'RelTol', 1e-20);
    while abs(cur_CDF-CDF_ary(i)) > min([CDF_ary, 1-CDF_ary, 10^-6])
        if cur_CDF > CDF_ary(i)
            max_X = cur_X;
        else
            min_X = cur_X;
        end
        cur_X = (min_X + max_X) / 2;
        cur_CDF = integral(@(y) pdf_func(y), -inf, cur_X, 'RelTol', 1e-20);
    end
    x_ary(i) = cur_X;
    
    
end

end