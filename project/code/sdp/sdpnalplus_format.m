function [blk, At, C, b, L, U, Bt, l, u, OPTIONS] = sdpnalplus_format(csp)
% Input arguments.
%   
%   csp - a "CSP" type object representing a well defined Constraint
%   Satisfaction Problem.
%
% Output arguments are for use in calls to SDPNAL+, in the following way:
%
% [obj,X,s,y,S,Z,y2,v,info,runhist] = ...
%  sdpnalplus(blk,At,C,b,L,U,Bt,l,u,OPTIONS);
%
% A limited description of some of these parameters is given below;
% see SDPNAL+ documentation for details.
%

Bt = []; % should never need to be changed.
l = []; % should never need to be changed.
u = []; % should never need to be changed.
OPTIONS.tol = 10^(-6); % default is 10^(-6)
OPTIONS.tolADM = 10^(-4); % default is 10^(-4)
OPTIONS.maxiter = 500; % 1000 is the default
OPTIONS.maxiterADM = 200; % the default is "near" 200.
OPTIONS.printlevel = 1; % 1 is the default
OPTIONS.stopoption = 0; % 1 implies "stop the solver if stagnation ocurs."; can set to 0.
OPTIONS.AAtsolve = 2; % 2 (the default) is to use a Cholesky decomposition; can set to 1.

% If the kth block X{k} of the variable X is a nonnegative vector block with
% dimension nk, then we set
% blk{k,1} = 'l', blk{k,2} = nk,
% At{k} = [nk  m sparse], Bt{k} = [nk  p sparse],
% C{k}, L{k}, U{k}, X{k}, S{k}, Z{k} = [nk  1 double or sparse].
%
% If the jth block X{j} of the variable X is a semidenite block consisting
% of a single block of size sj , then the content of the jth block is given
% as follows: blk{j,1} = 's', blk{j,2} = sj ,
% At{j} = [sjBAR  m sparse ], Bt{k} = [sjBAR  p sparse],

prob_vars = 1:(csp.numVariables);
dom = csp.domain; 
D = length(dom); 
N = D * csp.numVariables;
sig_space = sets2space({prob_vars, dom});

ALcell = cell(csp.numConstraints,1);
AScell = cell(csp.numConstraints,1);
Ccell = cell(csp.numConstraints,1);
lambda_dim = 0;
m = 0;
b = [];
tempC = [];
display(sprintf('Building individual-constraint matrices...'));
for i = 1:csp.numConstraints
    [ALcell{i}, AScell{i}, Ccell{i}] = ...
        constructConstraintsSDPMatrices(csp, i);
    if (mod(i-1,100) == 0)
        display(strcat('progress_',num2str(i/csp.numConstraints))); 
    end
    lambda_dim = lambda_dim + size(ALcell{i},2); % count number of columns
    nr = size(ALcell{i},1); % number of new rows
    m = m + size(ALcell{i},1); % total rows so far
    b = [b, 1, zeros(1,nr-1)]; % b won't be done constructing until later in this function.
    tempC = [tempC, Ccell{i}];
end
display(sprintf('Done building individual-constraint matrices. \n'));
C{1} = -tempC'; % we have a "maximization" objective, but SDPNAL+ assumes 
% a minimization objective (hence the negative sign).
C{2} = sparse(N,N); % the SDP block doesn't contribute to the objective.
C = C';

U{1} = sparse(ones(lambda_dim,1));
U{2} = sparse(ones(N));
U = U';
L{1} = sparse(zeros(lambda_dim,1));
L{2} = sparse(zeros(N));
L = L';

% do the SDP block first.
blk{2,1} = 's';
blk{2,2} = N;
s_bar = blk{2,2} * (blk{2,2} + 1)/2;

AS_force_zeros = {};
ct2 = 0;
for v = prob_vars
    for ell1 = 1:length(dom)
        for ell2 = (ell1+1):length(dom)
            val1 = dom(ell1);
            val2 = dom(ell2);
            % add a constriant forcing the corresponding SigmaVar to
            % zero.
            sig_row = tuple2linear(sig_space, [v, val1]);
            sig_col = tuple2linear(sig_space, [v, val2]);
            AS_force_zeros{ct2+1} = sparse([sig_row sig_col],...
                    [sig_col sig_row],...
                    [1 1], N, N);
            ct2 = ct2 + 1;
        end
    end
end
b = [b, zeros(1,ct2)];
m = m + ct2; % m == length(b)
b = b';

display(sprintf('Vectorizing main SDP constraint matrices...'));
AS = sparse(m , s_bar); % each matrix needs to be put in svec format.
count = 0;
for i = 1:csp.numConstraints
    if (mod(i-1,100) == 0)
        display(strcat('progress_',num2str(i/csp.numConstraints))); 
    end
    for z = 1:length(AScell{i})
        AS(count + 1,:) = svec(blk(2,:), AScell{i}{z});
        count = count + 1;
    end
end
display(sprintf('Done vectorizing main SDP constraint matrices. \n'));

display(sprintf('Vectorizing force-zero SDP constraint matrices...'));
for i = 1:ct2
    if (mod(i-1,100) == 0)
        display(strcat('progress_',num2str(i/csp.numConstraints))); 
    end
    AS(count + 1,:) = svec(blk(2,:), AS_force_zeros{i});
    count = count + 1;
end 
display(sprintf('Done vectorizing force-zero SDP constraint matrices. \n'));

% now do the linear block.
blk{1,1} = 'l';
blk{1,2} = lambda_dim;
AL = sparse(m, lambda_dim);
AL(1:size(ALcell{1},1), 1:size(ALcell{1},2)) = sparse(ALcell{1});
pr = size(ALcell{1},1); % previous row
pc = size(ALcell{1},2); % previous column
display(sprintf('Assembling linear constraint matrix...'));
for i = 2:csp.numConstraints
    if (mod(i-1,100) == 0)
        display(strcat('progress_',num2str(i/csp.numConstraints))); 
    end
    nr = size(ALcell{i},1); % number of new columns
    nc = size(ALcell{i},2); % number of new rows
    AL( (pr + 1):(pr + nr), (pc + 1):(pc + nc) ) = sparse(ALcell{i});
    pr = pr + nr;
    pc = pc + nc;
end
display(sprintf('Done assembling linear constraint matrix. \n'));

At{1} = AL';
At{2} = AS';
At = At';

end
