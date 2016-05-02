function sol = solve_csp_optimally(csp)

%%%%%%%%%%%%%%%%%%%$%%%%% define constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = csp.numVariables;
m = csp.numConstraints;
D = length(csp.domain);
N = n * D;

%%%%%%%%%%%%%% construct individual constraint matrices %%%%%%%%%%%%%%%

const_matr = cell(m,1);
C = cell(m,1);
for i = 1:m
    [const_matr{i}, C{i}] = constructConstraintsIPmatrices( csp, i );
end
m0 = 0; Ys = cell(m, 1); rhss = cell(m, 1);
for i = 1:m % this loop is just to reorganize data for easy access
    Ys{i} = const_matr{i}.Y_mat;
    rhss{i} = const_matr{i}.RHS;
    m0 = m0 + size(const_matr{i}.Y_mat,1);
end

%%%%%%%%%%%%%%%%%%%% combine constriant matrices %%%%%%%%%%%%%%%%%%%%%%%

behind_us = 0;
row_count_offset = 0;
row_ind = zeros(m0-m, 1);
col_ind = zeros(m0-m, 1);
allvals = zeros(m0-m, 1);
for i = 1:m
    nnzr = const_matr{i}.X_mat.num_rows ; % num nonzero rows
    nnzv = length(const_matr{i}.X_mat.colIdx); % num nonzero vals
    start_idx = behind_us + 1; % +1 for zero indexing
    end_idx = start_idx + (nnzv-1);
    tmp = imresize(2:(nnzr+1), [1 nnzv], 'nearest'); % "stretch" a vector.
    row_ind(start_idx:end_idx) = row_count_offset + tmp;
    col_ind(start_idx:end_idx) = const_matr{i}.X_mat.colIdx;
    allvals(start_idx:end_idx) = const_matr{i}.X_mat.vals;
    row_count_offset = row_count_offset + (nnzr + 1);
    behind_us = behind_us + nnzv;
end
Y = sparse(blkdiag(Ys{:}));
X = sparse(row_ind, col_ind, allvals, m0, N);
Xsum = kron(eye(csp.numVariables), ones(1,D));
Ysum = sparse(size(Xsum,1), size(Y,2));
full_mat = [Y, X; Ysum, Xsum];

%%%%%%%%%%%%%%%%%%%% define rest of Gurobi model %%%%%%%%%%%%%%%%%%%%%%%

f = [horzcat(C{:}), zeros(1,N)];
sense = char(1, size(full_mat,1));
equality_indices = find(ismember(X,zeros(1,N),'rows'))';
equality_indices = [equality_indices, (size(X,1)+1):(size(X,1)+D)];
sense(equality_indices) = '=';
sense(setdiff(1:(size(full_mat,1)),equality_indices)) = '<'; 
model.sense = sense;
model.modelsense = 'max';
model.A = full_mat;
model.rhs = [vertcat(rhss{:}); ones(csp.numVariables, 1)];
model.obj = f;
model.vtype = 'B';
params.method = -1;
save('about_to_solve_folded_csp.mat');

%%%%%%%%%%%%%%%%%%%% solve Gurobi model %%%%%%%%%%%%%%%%%%%%%%%

result = gurobi(model, params);
save('solved_folded_csp.mat')

%%%%%%%%%%%%%%%%%%%% record solution %%%%%%%%%%%%%%%%%%%%%%%

count = size(Y,2) +1; 
for i = 1:(csp.numVariables) 
    for j = 1:D 
        resMat(i,j) = result.x(count); 
        count = count + 1; 
    end
    assignment(i) = find(resMat(i,:));
end

sol.gurobi_output = result;
sol.assn = assignment;

end

