function [ muv, lambda, cvx_optval ] = lprelax( problem )

%% Simple LP Relaxation

%% Defining the CSP / recalling it

V = problem.numVariables; % Number of variables
D = length(problem.domain);% Size of domain, assuming <=10
arity = problem.arity;
constraints = problem.constraints; % Load the constraints
num_c = problem.numConstraints; % Number of constraints
M = D^arity; % Maximum number of local assignments for any constraint, uniform bound makes life easy


%% Evaluating the satisfiability matrix

RLS = zeros(num_c,M ); % RLS is consistent with R_i(L(S_i))
% Indicator for satisfiability of "i"th constraint by the "j"th local
% assignment. "M" number of assginments is an overkill but the time I spent
% couldn't figure out a smarter way.

ML = false(num_c, M, V, D); % ML stands for mu locations
% indexed by i, j, v, l 
% For each constraint "i", for each local assignment "j", for each variable
% "v", for each l=0,..., D-1; assign ML(i, j, v, l+1) = true if "j"th
% local assignment assigns variable "v", the value "l". Will come in handy
% while imposing the first constraint in the report for the LP relzaxation

% Now evaluate the satisfiability of all local assignments. That is
% evaluate RLS matrix (num_c constraints, M (max) possible local
% assignments). 

for i=1:num_c 
    constraint_i = constraints{i}; % loading the "i"th constraint
    scope_i = constraint_i.scope; % the variable indices for the constraint
    num_scope_i = length(scope_i); % cardinality of scope
    num_local_assgn_Ci = D^num_scope_i; % number of local assignments
    
    %Note that we will evaluate only n=num_local_assgn_Ci number of RLS
    %entries for the "i"th constraint. Rest "M-n" will be set to zero so
    %that they don't contribute to the objective and will be zero at
    %optimum as at optimal there is no point putting weights on constraints
    %with RLS entry zero. 
    for j=1:num_local_assgn_Ci
        
        x = zeros(V, 1); %Total Assignment(we can use simply local as well)
        s = dec2base(j-1, D, num_scope_i);
        % For D<=10, all digits are of the form 0, 1, ..., 9 but for D >=11
        % it starts using A, B, ..., Z and thus dec2base exists for D<=36
        % and is usable for the current logic for D<=10
        % The third argument of the function tells at least how many digits
        % should be present in the converted string.
        
        % s =[repmat('0', 1,  num_scope_i- length(s)), s];
        % The above two lines are in some sense a brute force to index the
        % various local assignments. First line finds a binary
        % representation for the current assignment and the second line
        % apends it with necessary number of zeros. Say for constraint with
        % scope x2 x3 x4, we convert 3 to binary we will get 11, so we
        % append a zero in front and get 011. This is a string so we
        % convert it to number term by term and assign to corresponding
        % index in x. 
        
        
        % That is if scope = [2 3 4], then s=011; and now we assign
        % x(scope(1)) = str2double(s(1)) and this is equivalent to
        % assigning x(2) = 0. 
        % Also, we need to set the corresponding indicator in ML as true.
        % Corresponding entry will have index as "i" for constraint, "j"
        % for the assignment, scope_i(k) for the variable and the
        % assignment value+1 for the last index.
        
        
        for k = 1:num_scope_i
            x(scope_i(k)) = str2double(s(k));
            ML(i, j, scope_i(k), x(scope_i(k))+1) = true;
            % Can write also as:
            % idx = sub2ind(size(ML), i, j, scope_i(k), x(scope_i(k))+1);
            % ML(idx)=true;
        end
        
        % Now evaluate the "i-j"the entry of RLS.
        RLS(i, j) = constraint_i.evaluate(x);
    end
end

W = diag(problem.weights); % weight for constraints num_c \times num_c

% The next few lines explain how we succinctly express our objective
% function using the varibles RLS, W and the LP variable Lambda that we
% will define in the next section.
%
% The objective can be expressed as 
% \sum_{i} weights_i \sum_{j} \lambda_{ij} RLS_{ij} which is same as 
% \sum_i \sum_j \lambda_{ij} (weights_i*RLS_{ij})

% Now if we define W = diag(weights), then for WRLS = W \times RLS, we have
% WRLS_{ij} = weights_i * RLS_{ij}.  And thus the objective becomes 
% \sum_i \sum_j \lambda_{ij} WRLS_{ij}.

% Tr(A'*B)  = \sum_{j} (A'*B)_{jj} = \sum_{j} \sum_{i} A'_{ji} B_{ij}
%           = \sum_{j} \sum_{i} A_{ij} B_{ij} 

% Thus our objective is same as Tr(\lambda' * WRLS)=Tr(\lambda' * W * RLS).

%% CVX Code for Optimization 

% muv(l) is a set of |V| |D| number of variables, we can represent
% them as a matrix with a row corresponding to a variable. So we will
% define a set of variables as a matrix of size V by D called muv (mu is
% reserved for MATLAB). The jth column will correspond to v = j-1.

% Lambda represents the dummy variable for all possible (M) local 
% assignments to the (num_c) constraints.

echo on

cvx_begin
    
    variable lambda(num_c, M)
    variable muv(V, D)
    maximize ( trace(lambda' * W * RLS ));
    subject to
        0<= lambda <=1;
        0<= muv <=1;
        for v = 1 : V
            sum(muv(v, :)) == 1;
        end
        for c_i=1:num_c
            sum(lambda(c_i, :)) == 1;
            for v = constraints{ c_i }.scope'
                for l=1:D
                    sum(lambda(c_i,  ML(c_i, :, v, l))) == muv(v, l);
                end
            end     
        end     
cvx_end

echo off
