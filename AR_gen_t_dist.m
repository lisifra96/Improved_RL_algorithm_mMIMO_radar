function c_norm=AR_gen_t_dist(N,Ntrans,rho,sigma2_c,lambda,scale)
%   This function generates an array of elements of an AR process with
%   t-distributed innovation process
%   INPUT PARAMETERS
%   N           = number of samples to generate
%   p           = aray containing the coefficients of the AR process
%   sigma2_c    = Variance of the generated samples
%   lambda      = parameter of the t-distribution
%   scale       = parameter of the t-distribution
%   OUTPUT PARAMETER
%   c           = array of elements of one realization of the process 

K = N+Ntrans;

% Innovations
w = sqrt(1/2)*(randn(K,1)+1j.*randn(K,1));
R = gamrnd(lambda,scale,K,1);
innov = w./sqrt(R);                             % T=Z/sqrt(GAMMA) This is a way to generate a t-distribution

cl = filter(1,rho,innov);                       % AR expression
c = cl(1+Ntrans:end);                           % Removal of the transient

sigma2_c_est = mean(abs(c).^2);                 % sample variance

c_norm = c*sqrt(sigma2_c/sigma2_c_est);         % scale to obtain the desired variance

end