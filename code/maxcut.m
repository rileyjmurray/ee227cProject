% Relation for the maxcut
r = @(x) x(1)~=x(2); 

% Creating the graph
n = 10; % number of vertices.
e = 20; % number of edges

% This rnd2 is simply doing - second = (second == first? n : second)
% It is to avoid any self edges in the graph that I generate randomly.
rnd2 = @(first, second, n)(second*(first~=second)+n*(first==second)); 
E = zeros(e,2);
for i=1:e
    first = randi(n); 
    second = randi(n-1); 
    second = rnd2(first, second, n);
    E(i,:) = [first, second];  
end

%{
E = (randi(n, e, 2)); % This way there are self pairs!
eZ = size(E); 
e = eZ(1); 
%}


%Weights, with min threshold and normalization
eps= 0.01; 
w = eps+rand(e, 1); 
w = w/sum(w); % This still gives an error in the error check in CSP.m

% Constraint Allocation in the format given
C = cell(e, 1); 
for i=1:e
    C{i} = Constraint(E(i,:), @(x) r(x)); 
end

% A random assignment - each variable taking value 0 or 1 independently
% with probability half.

assignment = (rand(n,1) >=0.5); 

maxCutProblem = CSP(w, C); 
maxCutProblem.evaluateObjective(assignment)
