function [new_vec,fixed_RN,fixed_RS,fixed_ZN,fixed_ZS] = fixpoles(old_vec,fixed_dim)
% Copyright @ Ruowen Liu, July 2023

if fixed_dim > 0
    design_vec_half_dim = length(old_vec)/2;
    fixed_RN = old_vec(1:fixed_dim);
    fixed_RS = old_vec(design_vec_half_dim-fixed_dim+1:design_vec_half_dim);
    fixed_ZN = old_vec(design_vec_half_dim+1:design_vec_half_dim+fixed_dim);
    fixed_ZS = old_vec(2*design_vec_half_dim-fixed_dim+1:2*design_vec_half_dim);
    new_vec = old_vec([(fixed_dim+1:design_vec_half_dim-fixed_dim)';(design_vec_half_dim+fixed_dim+1:2*design_vec_half_dim-fixed_dim)']);
else
    new_vec = old_vec;
    fixed_RN = []; fixed_RS = []; fixed_ZN = []; fixed_ZS = [];
end