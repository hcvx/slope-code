function absy_array = inverse_CDF_absy(p_array, list_dim, para)
size_p_array = size(p_array);
absy_array = zeros(size_p_array);


% sigma_H = para.sigma_H;
% sigma_X0 = para.sigma_X0;
% sigma_H1 = (sigma_H^2+sigma_X0^2)^0.5;
% mu_H1 = para.mu_X0;
% 
% list_dim_H0 = fix(list_dim/2);
% list_dim_H1 = list_dim - list_dim_H0;
% 
% absy_array_H0 = linspace(0, 5*sigma_H, list_dim_H0);
% absy_array_H1 = linspace(max(5*sigma_H, mu_H1-5*sigma_H1), mu_H1+5*sigma_H1, list_dim_H1+1);
% absy_array_H1 = absy_array_H1(2:list_dim_H1+1);
% 
% absy_array1 = [absy_array_H0 absy_array_H1];
% absy_array1 = sort(absy_array1, 'ascend');

pdf_y = @(y) pdf_of_X(y, para);
pdf_absy = @(y) pdf_y(y) + pdf_y(-y);
max_absy = invCDF_of_absX(1-10^-6, pdf_absy);
absy_array1 = linspace(0, max_absy, list_dim);




%absy_array1 = linspace(0, max_absy, list_dim);
CDF_list = CDF_of_absX(absy_array1, para);
CDF_list(list_dim) = 1;

%% new
% partition 
partition_length = 1024;
length_p_array = numel(p_array);
num_partition = ceil(length_p_array/partition_length);
for k=1:num_partition
    cur_partition_index = (k-1)*partition_length+1:min(length_p_array, k*partition_length);
    cur_p_array_partition = p_array(cur_partition_index);
    cmp_mtx = CDF_list'*ones(1,numel(cur_p_array_partition))<ones(list_dim,1)*cur_p_array_partition;
    absy_array1_index = min(length(absy_array1), sum(cmp_mtx,1)+1);
    absy_array(cur_partition_index) = absy_array1(absy_array1_index);
end
%% old
% for i=1:length(p_array)
%     absy_array(i) = absy_array1(min(length(absy_array1), sum(CDF_list<p_array(i))+1));
% end
end