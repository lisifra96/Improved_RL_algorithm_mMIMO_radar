function [obj, C] = Alg2v1(m,A,Pmax)
%% A is a matrix of M*N
%% Pmax is the maximum power
%% m is the dimension of the matrix C 
%% Intialize a feasible C
Cint = (rand(m) + 1i*rand(m));
for m1 = 1:m
    Cint(m1,:) = sqrt(Pmax/m)*Cint(m1,:)/norm(Cint(m1,:));
end
criteria = Inf;
%eps1 = 1e-10;
eps1 = 0.01;
tauC = 1e-6;
mul1 = 1.8;
C0 = Cint;
gamma1 =1;
cnt1 = 1;
N = size(A,2);
temp=zeros(N,1);
for n=1:N
    temp(n) = sum_square_abs(transpose(A(:,n))*C0);
end
obj(cnt1) = min(temp);
%% 

while(criteria > eps1)
    cvx_begin quiet                     % Prevents the model from producing any screen output while it is being solved.
    cvx_solver mosek
    
    % Variable declaration
    variable C(m,m) complex             % C is a mxm complex matrix
    variable t nonnegative              % t is a scalar non negative real number
    
    expression f
    maximise(t - tauC*norm(C-C0,'fro')) % Declaration of the function to maximise
    subject to 
    for n=1:N % N is the number of thetas
       derv1 = (2*conj(A(:,n))*((A(:,n).')*C0));
        -sum_square_abs((A(:,n).')*C0) - real(trace(((derv1)'*(C-C0)))) <= -t;
    end
    %% Power constraints
    % f = 0;
    % for m1=1:m
    % f = f+sum_square_abs(C(m1,:)); 
    % end
    % f<= Pmax;
    norm(C,'fro')<=sqrt(Pmax);
    cvx_end
    %% Update C0
    C0 = C0 + gamma1*(C-C0);
    gamma1 = gamma1*(1-1e-3*gamma1);
    cnt1 = cnt1+1; 
    obj(cnt1) = t;
    criteria = abs(obj(cnt1)-obj(cnt1-1));
    if(tauC < .5)
        tauC = tauC*mul1;
    end
    if(cnt1>50)
        break;
    end
end
end