function [ Alambda, Asigma, C ] = constructConstraintsSDPMatrices( csp, constr )
%
%   csp - a CSP object for the problem we're constructing an SDP for.
%
%   constr - the number of the constraint whose matrices Alambda and Asigma
%   we are trying to build.
%
%   Alambda - a matrix so that Alambda * lambda produces the
%   left-hand-sides of constraints (4.2) and (4.3) in Basic SDP (for
%   constraint "constr"). (Where "lambda" is a vector where the j^th entry
%   corresponds to the j^th local assignment for constraint "constr". I.e.,
%   "lambda" is of length "number of possible local assignments for this
%   constraint".
%
%   Asigma - a cell array of matrices so that Asigma{i}*sigmaVar is the
%   right-hand-side of constraint (4.3).
%       
%       In Basic SDP, sigmaVar has 4 indicies (e.g. sigmaVar(v,l,v',l')).
%       Within the solver, sigmaVar will only have 2 indicies (e.g.
%       sigmaVar(a,b), where "a" uniquely identifies (v,l) and "b" uniquely
%       identifies (v',l')). We acheive this by indexing sigmaVar as
%       sigmaVar((v * |D| + l),(v' * |D| + l')) (i.e. for a fixed variable,
%       entries corresponding to different values of that variable are 
%       consecutive; for a fixed value, entries correspond to different
%       variables are strided).
%   


% get "global" information about the larger problem
prob_vars = 1:(csp.numVariables);
dom = csp.domain; D = length(dom); N = D * csp.numVariables;
sig_space = sets2space({prob_vars, dom});

% get the information associated with this constraint
constr_obj = csp.constraints{constr};
scope = constr_obj.scope;
arity = length(scope);
lambda_space = sets2space({dom},arity);

% Alambda will have a number of columns equal to the number of local
% assignents for constraint "constr". The number of rows is harder to
% specify ahead of time.
num_cols = D^arity;
Alambda(1,:) = ones(1,num_cols);
temp = zeros(1,2);
Asigma{1} = sparse(N, N);
count = 2;
for a = 1:arity
    for b = a:arity % used to be 1:arity; think maybe should have that some times and not others? Shouldn't have for Ramsey???
        var1 = scope(a);
        var2 = scope(b);
        localassn = NaN(1,arity);
        if (var1 == var2)
            % define a row for all ell \in D.
            for ell = dom
                % update Asigma
                sig_row = tuple2linear(sig_space, [var1, ell]);
                sig_col = sig_row;
                Asigma{count} = sparse(sig_row, sig_col, -1, N, N);
                % temp(count,:) = [min([sig_row, sig_col]), max([sig_row, sig_col])]; for debugging.
                % update Alambda
                localassn(a) = ell; % localassn(b) = ell; is redundant
                Alambda(count,:) = zeros(1,num_cols);
                Alambda(count,inctuple2alllinear(lambda_space, localassn)) = 1;
                
                count = count + 1;
            end
        else
            % define a row for all (l,l') \in D \cross D
            for l1 = dom
                for l2 = dom
                    % update Asigma
                    sig_row = tuple2linear(sig_space, [var1, l1]);
                    sig_col = tuple2linear(sig_space, [var2, l2]);
                    Asigma{count} = sparse([sig_row sig_col],...
                        [sig_col sig_row],...
                        [-1 -1], N, N);
                    % temp(count,:) = [min([sig_row, sig_col]), max([sig_row, sig_col])]; for debugging.
                    % update Alambda
                    localassn(a) = l1;
                    localassn(b) = l2;
                    Alambda(count,:) = zeros(1,num_cols);
                    Alambda(count,inctuple2alllinear(lambda_space, localassn)) = 1;
                    count = count + 1;
                end
            end
        end
    end
end
Asigma = Asigma'; 

% C will be a vector so that C * lambda = sum_i sum_L w_i * R_i[L] * lambda_i[L] 
C =  zeros(num_cols, 1); % is a column vector with number of rows that matches "lambda"
for i = 1:num_cols
    tuple = linear2tuple(lambda_space, i);
    C(i) = csp.constraints{constr}.evaluate227c(tuple);
end
C = csp.weights(constr) * C';

end

