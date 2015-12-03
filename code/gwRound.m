%% gwRound(CSP)
    %   SDP produces sigmaVar and lambda
    %        sigmaVar : rows / cols are  indexed like 
    %        (v_1,l_1),(v_1,l_2),...,(v_1,l_D),(v_2,l_1)....
    % returns an assignment of variables for "CSP"
    
%% intialize and solve a particular example
worldColoring;
problem = africaProblem;
[sigmaVar, lambda] = constructAndSolveSDP(problem);
U = chol(sigmaVar);

%% Rounding
n = size(U,1); % number of rows
% generate normal vector for random hyperplane in R^n
v = 2 * (rand(1,n) - 0.5 * ones(1,n));
v = v / norm(v);
assignment = zeros(1,problem.numVariables);
for i = 1:problem.numVariables
    assignment(i) = sign(U(:,i) * n);
end