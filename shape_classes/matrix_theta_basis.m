function THETA_bas = matrix_theta_basis(t, NL, L)
% This computes the matrix THETA_bas such that
% [tht(t)_1, ..., tht(t)_{N-2}] = THETA_bas * design_vec_bas
% Here, N = NL+5 = length(design_vec)/2 + 2
% dimensions:
% t                    is (p*np) by 1
% THETA_bas            is (p*np) by 2*(N-2)
% THETA_bas(t)_j       is (p*np) by 1 (j-th column)

matR = [-26, -66, -26, -1, zeros(1, NL-1); ...
    eye(NL+3); ...
    zeros(1, NL-1), -1, -26, -66, -26];
matZ = 1i*[-10, 0, 10, 1, zeros(1, NL-1); ...
    eye(NL+3); ...
    zeros(1, NL-1), 1, 10, 0, -10];
THETA_bas = repmat( B0_original( (t*NL/L) - ((-5):(NL-1)) ), 1, 2 )* ( blkdiag(matR,matZ) );
end
