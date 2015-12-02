%% gwRound(CSP)
    % assumes a binary 2-CSP has been appropriately instantiated.
    %   this instance is simply "CSP"
    % assumes the SDP for this CSP has been solved.
    %   the last line sdprelax has the cholesky decomposition of sigmaVar
    %   this is stored in a variable "U".
    %   U = chol( sigmaVar );
    % returns an assignment of variables for "CSP"
    
%% Rounding
n = CSP.numVariables;
% generate random hyperplane
v = 2 * (rand(1,n) - 0.5 * ones(1,n));