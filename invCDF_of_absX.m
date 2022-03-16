function x_ary = invCDF_of_absX(CDF_ary, pdf_absX_func)
% find the quantile for a positive distribution
% pdf_absX_func: pdf of the positive distribution
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
    
    min_X = 0;
    max_X = 2;
    CDF_max_X = integral(@(y) pdf_absX_func(y), 0, max_X, 'RelTol', 1e-10);
    while CDF_max_X < max([CDF_ary, 1 - 10^-4])
        max_X = max_X * 1.5;
        CDF_max_X = integral(@(y) pdf_absX_func(y), 0, max_X, 'RelTol', 1e-10);
    end
    
    cur_X = (min_X + max_X) / 2;
    cur_CDF = integral(@(y) pdf_absX_func(y), 0, cur_X, 'RelTol', 1e-10);
    while abs(cur_CDF-CDF_ary(i)) > min([CDF_ary, 1-CDF_ary, 10^-4])
        if cur_CDF > CDF_ary(i)
            max_X = cur_X;
        else
            min_X = cur_X;
        end
        cur_X = (min_X + max_X) / 2;
        cur_CDF = integral(@(y) pdf_absX_func(y), 0, cur_X, 'RelTol', 1e-10);
    end
    x_ary(i) = cur_X;
end

end