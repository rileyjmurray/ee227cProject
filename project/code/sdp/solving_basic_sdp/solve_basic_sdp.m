function sol = solve_basic_sdp( varargin )
% The first argument must be a CSP object.
%
% All subsequent arguments come in name-value pairs, and currently only
% serve to specify options for SDPNAL+. In the future, name-value pairs can
% could be made more all-encompassing (particularly, we will allow the user
% to select different solvers).
%
if (isa(varargin{1},'CSP'))
    csp = varargin{1};
else
    error('First argument must be a CSP object.');
end

% "OPTIONS" returned by the function below uses default values.
[blk, At, C, b, L, U, Bt, l, u, OPTIONS] = sdpnalplus_format(csp);

for i = 2:2:(length(varargin)-1)
    switch varargin{i}
        case 'tol'
            OPTIONS.tol = varargin{i+1};
        case 'tolADM'
            OPTIONS.tolADM = varargin{i+1};
        case 'maxiter'
            OPTIONS.maxiter = varargin{i+1};
        case 'maxiterADM'
            OPTIONS.maxiterADM = varargin{i+1};
        case 'stopoption'
            OPTIONS.stopoption = varargin{i+1};
        case 'AAtsolve'
            OPTIONS.AAtsolve = varargin{i+1};
    end
end

id = datestr(datetime('now'),30);
tmp = computer;
save(strcat(tmp,id,'.mat')); % don't know exactly where this is going...

% [obj,X,s,y,S,Z,y2,v,info,runhist]  note which arguments are dropped
try
[obj,X,~,~,~,~,~,~,info,runhist] = ...
    sdpnalplus(blk,At,C,b,L,U,Bt,l,u,OPTIONS);
catch ME
    if(strcmp(ME.message(1:10),'mexbwsolve') || strcmp(ME.message(1:10),'mexfwsolve'))
        str = strcat({'SDPNAL+ exited with an error (likely due to'},...
            {' '},{'linearly dependent constraints). \n'},{' '},...
            {'Executing a subroutine to remove linearly dependent constraints'},...
            {' '},{'(this may take a while).'});
        warning(char(str));
        [AL2, AS2, b2] = remove_ld_rows(At{1}', At{2}', b);
        b = b2;
        AS = AS2;
        AL = AL2;
        At{1} = AL';
        At{2} = AS';
        [obj,X,~,~,~,~,~,~,info,runhist] = ...
            sdpnalplus(blk,At,C,b,L,U,Bt,l,u,OPTIONS);
    end
end

[L,D] = ldl(X{2});
M = L*sqrt(D);

sol.lambdaVar = X{1};
sol.sigmaVar = X{2};
sol.sigmaFactor = M;
sol.obj = obj;
sol.termcode = info.termcode;
sol.all_info = info;
sol.run_hist = runhist;

save(strcat(tmp,id,'.mat'));

end

