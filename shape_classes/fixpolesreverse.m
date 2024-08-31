function old_vec = fixpolesreverse(new_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS)
% Copyright @ Ruowen Liu, July 2023

fixed_dim = length(fixed_RN);
if fixed_dim == 0
    old_vec = new_vec;
else
    new_half_dim = length(new_vec)/2;
    old_vec = [fixed_RN; new_vec(1:new_half_dim); fixed_RS; fixed_ZN; new_vec(new_half_dim+1:end); fixed_ZS];
end
