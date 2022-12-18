function C = Closed_Form_W(A,Pmax)
%% A is a matrix of M*N
%% Pmax is the maximum power
%% m is the dimension of the matrix C 


[m, num_targets] = size(A);
G=zeros(m,m);

  for i=1:num_targets
    G=G+(conj(A(:,i))*(A(:,i)).');
  end

 tmp=(10^-10).*eye(m,m);
 G=G+tmp;
 G = (G+G')/2; % Be sure that the matrix is Hermitian
 G=(Pmax/trace(G))*G;
 C=sqrtm(G);

end