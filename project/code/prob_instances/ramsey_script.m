
%% Ramsey number instance
m = 3; % i.e., suppose we are looking for RAMSEY(3)
n = 13; % Is there a 2-coloring of K_10 such that there are no homogeneous sets of size 3?
domain = [0,1];
%% associate a variable index to each edge in K_n
edge_space = sets2space({1:n, 1:n}); % use special function for this space which skips self-edges
all_edges = nchoosek(1:n,2);
edge_linear2tuple = containers.Map('KeyType','double','ValueType','any');
for i = 1:length(all_edges)
   edge_linear2tuple(i) = all_edges(i,:);
end

%% define the set of constraints
all_mples = nchoosek(1:n,m);
relation = @(x) max(x)*(1-min(x)); % at least one 1 and at least one 0.
num_constrs = length(all_mples);
ramsey_constraints = cell(num_constrs, 1);
for i = 1:length(all_mples)
   scope = this_mples_edges(edge_space, all_mples(i,:));
   ramsey_constraints{i} = Constraint(scope, relation);
end
%% define the CSP
weights = 1/num_constrs * ones(num_constrs,1);
ramsey_csp = CSP(weights, ramsey_constraints, domain, n*(n-1)/2 );
%% solve the SDP
%sol = solve_basic_sdp(ramsey_csp, 'tol', 10^(-10));


    
    
    
    
    
    