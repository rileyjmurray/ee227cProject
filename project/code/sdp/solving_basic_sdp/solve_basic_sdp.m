function sol = solve_basic_sdp( varargin )
%
% ***************************** usage *************************************
%
% This function constructs Basic SDP for an input CSP in SDPNAL+ format,
% and then solves that SDP with SDPNAL+. Use this function for small CSP's,
% or when you do not have access to Matlab's Parallel Computing Toolbox.
%
% ************************** signature ************************************
%
%   sol = solve_basic_sdp_parallel_construct(csp), or
%
%   sol = solve_basic_sdp_parallel_construct(csp, Name1, val1, Name2, val2, ...)
%       "Name" is a setting for SDPNAL+ and "val" is the corresponding
%       value that the user would like that setting set to.
%
%   Name-value pairs are to specify settings for SDPNAL+.
%   Possible name-value pairs (name, example_value):
%          ('maxiter', 1000),
%          ('maxiterADM',500),
%          ('tol', 10^(-8)),
%          ('tolADM', 10^(-4)),
%          ('stopoption', 0),
%
%       ... see SDPNAL+ documentation for what each of these controls.
%
% **************************** inputs ************************************
%   
%   csp - a well defined CSP object.
%
% **************************** outputs ************************************
%
%   sol - a struct containing the following fields:
%
%       sol.lambdaVar - the linear variables in Basic SDP.
%       sol.sigmaVar - the matrix variable in Basic SDP.
%       sol.sigmaFactor - the Cholesky (or LDL') factorization of the SDP
%           matrix variable
%       sol.obj - the value of the SDP primal objective (obj(1)) and the
%           SDP dual objective (obj(2))
%       sol.termcode - the termination code of SDPNAL+.
%       sol.all_info - the "info" struct returned by SDPNAL+.
%       sol.run_hist - the "runhist" struct returned by SDPNAL+.
%   
%   NOTE : This function writes data to a file TWICE. 
%   
%       The first time it writes data to a file is after completing a call
%       to sdpnalplus_format. It saves the data at that point because calls 
%       to sdpnalplus_format may take a long time, and if SDPNAL+ exits 
%       with an error the user will likely want to refer back to the model
%       data without having to re-run sdpnalplus_format.
%
%       The second time it writes data to a file is after SDPNAL+ solves
%       the SDP. It saves data at this point because running SDPNAL+ will
%       take a long time for large problems, and we'd rather not lose
%       the solution because the user didn't store a variable properly
%       (e.g., sol = solve_basic_sdp_parallel_construct( csps{i} ), when
%       the user intended 
%       sols{i} = solve_basic_sdp_parallel_construct( csps{i} ) ).
%
% *************************************************************************

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

