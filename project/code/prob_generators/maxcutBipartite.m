function [ maxCutProblem ] = maxcutBipartite( n, e )
% generate random bipartite graph for max cut problem 
% n = number of vertices
% e = number of edges

%rng( seed )

% Relation for the maxcut
r = @(x) x(1)~=x(2); 

V = 1:n; 

V1 = unique(randi(n, round(n/2), 1)); 

V2 = setdiff(V, V1)';

E1 = randsample(V1, e, 1); 

E2 = randsample(V2, e, 1); 

E = [E1 E2]; 
E = unique(E,'rows');
e = size(E,1);

%Weights, with min threshold and normalization
eps= 0.01; 
w = eps+rand(e, 1); 
w = w/sum(w); % This still gives an error in the error check in CSP.m

% Constraint Allocation in the format given
C = cell(e, 1); 
for i=1:e
    C{i} = Constraint(E(i,:), @(x) r(x)); 
end

maxCutProblem = CSP(w, C, [0 1], n); 


