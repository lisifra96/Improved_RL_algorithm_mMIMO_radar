 function W_T = Wald_test(N,v,x,hat_alpha,l,n2)
%   This function computes the Wald Test decision statistic Λ as described
%   in "Massive MIMO Radar for Target Detection" sec.IIIB.
%   INPUT PARAMETERS
%   N           = number of elements in v and x
%   v           = view expression (12) in the paper
%   x           = view expression (12) in the paper
%   hat_alpha   = view expression (15) in the paper
%   l           = view expression (20) in the paper
%   n2          = v'*v
%   OUTPUT PARAMETERS
%   W_T         = Wald test decision statistics

B_est_mod = B_matrix_est(N,v,x,hat_alpha,l);
    
% Wald Test

W_T = 2 * n2^2 * (abs(hat_alpha)^2) /B_est_mod;

end

