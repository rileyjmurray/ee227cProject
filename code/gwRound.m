function [assignment] = gwRound(sigmaVar)
%% gwRound(CSP)
    %   SDP produces sigmaVar and lambda
    %        sigmaVar : rows / cols are  indexed like 
    %        (v_1,l_1),(v_1,l_2),...,(v_1,l_D),(v_2,l_1)....
    % returns an assignment of variables for "CSP"
    
%% intialize and solve a particular example

U = chol(sigmaVar);

%% Rounding
n = size(U,1); % number of rows

% generate normal vector for random hyperplane in R^n
mu = zeros(n,1);
covar = eye(n); 
v = mvnrnd(mu,covar,1); 
% v = 2 * (rand(1,n) - 0.5 * ones(1,n));

assignment = zeros(1, n);
for i = 1:n
    assignment(i) = sign(v*U(:,i)) ;
end
