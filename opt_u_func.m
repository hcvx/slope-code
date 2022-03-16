function output_u_array = opt_u_func(absy_array, u_array, input_absy_array)
% opt_u_func: for the proximal function prox() described by [absy_array, u_array]
% compute output_u_array = prox(input_absy_array)

size_output = size(input_absy_array);
output_u_array = zeros(size_output);

for i=1:length(input_absy_array)
    
    cur_idx = max(1,sum(absy_array<input_absy_array(i)));
    if cur_idx<length(absy_array)
        r = (absy_array(cur_idx+1)-input_absy_array(i))/(absy_array(cur_idx+1)-absy_array(cur_idx));
        output_u_array(i) = r*u_array(cur_idx)+(1-r)*u_array(cur_idx+1); 
    else
        output_u_array(i) = input_absy_array(i) - (-u_array(length(u_array))+absy_array(length(u_array)));
    end
       
end

end