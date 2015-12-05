function [sigmaVar, cvx_optval]= gwSDP(problem)

n = problem.numVariables;
C = problem.constraints;
num_c = length(C); 
w = problem.weights;

W = zeros(n, n); 
for i=1:num_c
    constraint_i = C{i}; % loading the "i"th constraint
    scope_i = constraint_i.scope; % the variable indices for the constraint
    W(scope_i(1), scope_i(2)) = W(scope_i(1), scope_i(2))+w(i); 
    W(scope_i(2), scope_i(1)) = W(scope_i(2), scope_i(1))+w(i); 
end

echo on 

cvx_begin
    variable sigmaVar(n,n) symmetric semidefinite
    maximize (0.5*(sum(w)-trace(W'*sigmaVar))); 
    subject to 
        diag(sigmaVar) == 1;
cvx_end

echo off


