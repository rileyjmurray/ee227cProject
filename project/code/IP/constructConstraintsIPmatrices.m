function [ const_matr, C] = constructConstraintsIPmatrices( csp, constr )
%
%   csp - a CSP object for the problem we're constructing an IP for.
%

% get "global" information about the larger problem
dom = csp.domain; D = length(dom); N = D * csp.numVariables;

% get the information associated with this constraint
constr_obj = csp.constraints{constr};
scope = constr_obj.scope;
arity = length(scope);
lambda_space = sets2space({dom},arity);

num_loc_assn = D^arity;
const_matr.Y_mat = sparse([ones(1,num_loc_assn); eye(num_loc_assn)]); % done.
const_matr.RHS = zeros(num_loc_assn + 1,1); const_matr.RHS(1) = 1; % done.

const_matr.X_mat.num_rows = num_loc_assn;
const_matr.X_mat.colIdx = zeros(arity * num_loc_assn,1);
const_matr.X_mat.vals = (-1/(arity-0.5))*ones(arity * num_loc_assn,1);
% const_matr.fullX = zeros(num_loc_assn, N); debugging

count = 1;
C =  zeros(num_loc_assn, 1); % is a column vector with number of rows that matches "lambda"
for i = 1:num_loc_assn
    tuple = linear2tuple(lambda_space, i);
    locs = zeros(1,arity);
    for j = 1:arity
       locs(j) = (scope(j)-1)*D + find(dom == tuple(j));
       const_matr.X_mat.colIdx(count) = locs(j);
       count = count + 1;
    end
    % const_matr.fullX(i,locs) = (-1/(arity-0.5)); debugging
    C(i) = csp.constraints{constr}.evaluate227c(tuple);
end
C = csp.weights(constr) * C';
end

