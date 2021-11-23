function [B_est] = B_matrix_est(N,v,x,alpha,l)
%   This function computes the B matrix that you can find in the article
%   "Massive MIMO Radar for Target Detection" sec.IIIA. B is N^(-1)*B_est in the paper.
%   The signal model in the presence of target is x=alpha*v+c.
% INPUT VARIABLES
% N:        Integer
% v:        array with N elements
% x:        array with N elements
% alpha:    scalar
% l:        truncation lag (l=ceil(N^(1/4)));
% OUTPUT VARIABLES
% B_est:    estimate of B (scalar)

hat_c = x - alpha*v;
y=hat_c.*conj(v);
yc=conj(y);

B_est = sum(abs(y).^2);
for m = 1:l
    B_est = B_est + 2*sum(real(y(m+1:N) .* yc(1:N-m)));
end

end
