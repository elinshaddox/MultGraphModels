% This script generates 4 random matrices similar to those used 
% in simulation 1 i.e. with some overlapping and some differential
% structure

% Create true precision matrices A1, A2, A3, and A4

% K is the number of subgroups
K = 4;

% p is the number of variables
p = 20;

% A1 is an AR(2) graph
A1 = toeplitz([1, 0.5, 0.4, zeros(1, p - 3)]);
n_edges = (sum(sum(A1 ~= 0)) - p) / 2;
n_possible = p * (p - 1) / 2;

Omega1 = toeplitz([1, 0.5, 0.4, zeros(1, p - 3)]);
tedge = p*(p-1)/2; beta = 0.5; 
indmx = reshape([1:p^2],p,p); 
upperind = indmx(triu(indmx,1)>0); 

% Locations of all nonzero entries of A1 above the diagonal
[rposA1, cposA1] = find(triu(A1) - eye(p));

% Locations of all zero entries of A1 above the diagonal
[rzeroA1, czeroA1] = find(triu(A1 == 0));

 %Randomly generating Omega2 with 70% of nodes of Omega1
 % i.e. randomly changing a fixed fraction of nonzero edges
 
 % desired fraction to set to zero
 f = 0.3; 
 
 %find nonzero edges in Omega1
 nonzero = find(Omega1);
 n = round(f*length(nonzero)); 

nonzero_change = randsample(nonzero,n);
Omega2=Omega1;
% apply changes in those entries
Omega2(nonzero_change) = 0; 

% To make positive definite
% If you don't get a positive-definite matrix, you can play with
% the parameters c and d to get one.
Network1=Omega2;
c = 3.4; d = 0.6;
AA = Network1.*c;
AA = (AA+AA')>0;
AA = AA-diag(diag(AA))+eye(p)*d;
AA =  AA./(ones(p,1)*(1.4*sum(abs(AA))))';
AA = (AA + AA')./2;
AA = AA-diag(diag(AA))+eye(p)*d;
Omega2 = AA;

% To check matrix is positive definite
% If the input matrix is not positive definite, then "p" will be a positive integer
[~,pos]=chol(AA)

% To generate Omega 3 by changing five nonzero entries of Omega 1 to zeroes
nonzero = find(Omega1);
zero = find(~Omega1);
nonzero_change = randsample(nonzero,5);
zero_change = randsample(zero, 5);
Omega3=Omega1;
Omega3(nonzero_change) = 0; 
Omega3(zero_change) = .6;
Network2=Omega3;
c = 3.4; d = 0.6;
AA = Network2.*c;
AA = (AA+AA')>0;
AA = AA-diag(diag(AA))+eye(p)*d;
AA =  AA./(ones(p,1)*(1.4*sum(abs(AA))))';
AA = (AA + AA')./2;
AA = AA-diag(diag(AA))+eye(p)*d;
Omega3 = AA;

% To check that Omega3 is positive definite
[~,pos]=chol(Omega3)

% To generate Omega4
b_prior = 3; D_prior = eye(p); n = 3*p; b_post = b_prior+n;
A = toeplitz([1,0.5,zeros(1,p-2)]); 
A(1,p)=0.4; A(p,1) = 0.4;
Omega4 = A;
[~,pos]=chol(Omega4)




% Save out these matrices to use as input to simulation script
csvwrite('A1_sim1.csv', Omega1)
csvwrite('A2_sim1.csv', Omega2)
csvwrite('A3_sim1.csv', Omega3)
csvwrite('A4_sim1.csv', Omega4)
