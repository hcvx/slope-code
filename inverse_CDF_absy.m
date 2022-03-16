function absy_array = inverse_CDF_absy(p_array, list_dim, para)
size_p_array = size(p_array);
absy_array = zeros(size_p_array);


pdf_y = @(y) pdf_of_X(y, para);
pdf_absy = @(y) pdf_y(y) + pdf_y(-y);
max_absy = invCDF_of_absX(1-10^-6, pdf_absy);
absy_array1 = linspace(0, max_absy, list_dim);

CDF_list = CDF_of_absX(absy_array1, para);
CDF_list(list_dim) = 1;


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

end