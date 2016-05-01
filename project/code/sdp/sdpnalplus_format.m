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
% A limited description of some of these parameters is given at the END of 
% this function; see SDPNAL+ documentation for details.
%

%%%%%%%%%%%%%%%%%%%%%%%%%% Set default options %%%%%%%%%%%%%%%%%%%%%%%%%%

OPTIONS.tol = 10^(-6); % default is 10^(-6)
OPTIONS.tolADM = 10^(-4); % default is 10^(-4)
OPTIONS.maxiter = 500; % 1000 is the default
OPTIONS.maxiterADM = 200; % the default is "near" 200.
OPTIONS.printlevel = 1; % 1 is the default
OPTIONS.stopoption = 0; % 1 implies "stop the solver if stagnation ocurs."; can set to 0.
OPTIONS.AAtsolve = 2; % 2 (the default) is to use a Cholesky decomposition; can set to 1.

%%%%%%%%%%%%%%%%% Give short names to important things %%%%%%%%%%%%%%%%%%

prob_vars = 1:(csp.numVariables);
dom = csp.domain; 
D = length(dom); 
N = D * csp.numVariables;
Nbar = N * (N + 1)/2;
sig_space = sets2space({prob_vars, dom});

ALcell = cell(csp.numConstraints,1);
AScell = cell(csp.numConstraints,1);
Ccell = cell(csp.numConstraints,1);
m = 0; bIdx = [];

%%%%%%%%%%%%%%% Get matrices for each "constraint" in CSP %%%%%%%%%%%%%%%%

display(sprintf('\n Building individual-constraint matrices...'));
for i = 1:csp.numConstraints
    [ALcell{i}, AScell{i}, Ccell{i}] = ...
        constructConstraintsSDPMatricesNoSvec(csp, i);
    if (mod(i-1,100) == 0)
         display(char(strcat({'progress'},{' '},{num2str(i/csp.numConstraints)}))); 
    end
    bIdx = [bIdx, m + 1];
    m = m + size(ALcell{i},1); % total rows so far
end
b = full(sparse(1, bIdx, 1, 1, m));
%%%%%%%%%%%%%%%% Construct force-zero constraint matrices %%%%%%%%%%%%%%%%

ASFZvalues = []; ASFZcolindices = []; lfz = 0;
for v = prob_vars
    for ell1 = 1:length(dom)
        for ell2 = (ell1+1):length(dom)
            val1 = dom(ell1);
            val2 = dom(ell2);
            % add a constriant forcing the corresponding SigmaVar to
            % zero.
            sig_row = tuple2linear(sig_space, [v, val1]);
            sig_col = tuple2linear(sig_space, [v, val2]);
            si = rowcol2svecidx(sig_row, sig_col);
            ASFZcolindices = [ASFZcolindices, si];
            if (si == 1 || si == Nbar)
                ASFZvalues = [ASFZvalues, 1];
            else
                ASFZvalues = [ASFZvalues, sqrt(2)];
            end
            lfz = lfz + 1;
        end
    end
end
ASFZ = sparse(1:lfz,ASFZcolindices,ASFZvalues,lfz,Nbar);
b = [b, zeros(1,lfz)];
% m = m + lfz; % m == length(b), don't need "m" anymore.

%%%%%%%%%%%%%%%% Assemble SDP variable constraint matrix %%%%%%%%%%%%%%%%

% display(sprintf('\n Assembling main SDP constraint matrices...'));
% AS = sparse(m , Nbar);
% count = 0; 
% for i = 1:csp.numConstraints
%     if (mod(i-1,100) == 0)
%          display(char(strcat({'progress'},{' '},{num2str(i/csp.numConstraints)}))); 
%     end
%     step = size(AScell{i},1);
%     AS((count + 1):(count + step),:) = AScell{i};
%     count = count + step;
% end
% step = lfz;
% AS((count + 1):(count + step),:) = ASFZ;
temp = [AScell; {ASFZ}];
AS = vertcat(temp{:});

%%%%%%%%%%%%%%% Assemble linear variable constraint matrix %%%%%%%%%%%%%%%

% AL = sparse(m, lambda_dim);
% AL(1:size(ALcell{1},1), 1:size(ALcell{1},2)) = sparse(ALcell{1});
% pr = size(ALcell{1},1); % previous row
% pc = size(ALcell{1},2); % previous column
% display(sprintf('\n Assembling linear constraint matrix...'));
% for i = 2:csp.numConstraints
%     if (mod(i-1,100) == 0)
%         display(char(strcat({'progress'},{' '},{num2str(i/csp.numConstraints)}))); 
%     end
%     nr = size(ALcell{i},1); % number of new columns
%     nc = size(ALcell{i},2); % number of new rows
%     AL( (pr + 1):(pr + nr), (pc + 1):(pc + nc) ) = sparse(ALcell{i});
%     pr = pr + nr;
%     pc = pc + nc;
% end
temp = blkdiag(ALcell{:});
lambda_dim = size(temp, 2);
AL = vertcat(temp,sparse(lfz,lambda_dim));

%%%%%%%%%%%%%%%%%%%%%%%%%% Return values %%%%%%%%%%%%%%%%%%%%%%%%%%

blk{2,1} = 's';
blk{2,2} = N;
blk{1,1} = 'l';
blk{1,2} = lambda_dim;
Bt = []; l = []; u = [];  % these should never be changed.
C{1} = -(horzcat(Ccell{:}))'; % we have a "max" objective, but SDPNAL+ assumes a "min".
C{2} = sparse(N,N); % the SDP block doesn't contribute to the objective.
C = C'; % This is only reshaping the cell array for SDPNAL+!
b = b'; % b should be a column vector.
At{1} = AL';
At{2} = AS';
At = At'; % This is only reshaping the cell array for SDPNAL+!
U{1} = sparse(ones(lambda_dim,1));
U{2} = sparse(ones(N));
U = U'; % This is only reshaping the cell array for SDPNAL+!
L{1} = sparse(zeros(lambda_dim,1));
L{2} = sparse(zeros(N));
L = L'; % This is only reshaping the cell array for SDPNAL+!
display(sprintf('\n Done converting to SDPNAL+ format. \n'));

end

% FROM THE SDPNAL+ guide:
%
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

