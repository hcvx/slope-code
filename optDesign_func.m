function [absy_array, u_array, sigma_H, tau] = optDesign_func(num_of_grid, para)

flag = 0;

sigma = para.sigma;

delta = para.delta;

iter = 0;
max_Iter = 100;
length_threshold = 10^-6;

sigma_z_max = (sigma^2 + para.sgl_sq/delta)^0.5;

cur_sigma_z_left = sigma;
cur_sigma_z_right = sigma_z_max;

while flag==0
    
    sigma_H = (cur_sigma_z_left + cur_sigma_z_right)/2;
    para.sigma_H = sigma_H;   

    [absy_array, u_array, f_eval, Expe_deri] = optProx_func(num_of_grid, para);
    
    prox_MSE = f_eval + para.sgl_sq;
    
    expected_MSE = (sigma_H^2-sigma^2)*delta;
    minimum_MSE = prox_MSE;
    
    if (expected_MSE<minimum_MSE)
        cur_sigma_z_left = sigma_H;
    else
        cur_sigma_z_right = sigma_H;
    end
    iter = iter + 1;
    
    if abs(cur_sigma_z_right-cur_sigma_z_left)<length_threshold || iter>max_Iter
        flag=1;
    end
end

aver_slope = Expe_deri;

tau = 1/(1-aver_slope/delta);
end
