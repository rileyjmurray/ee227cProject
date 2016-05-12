function sol = vfm(csp, epsilon)
%
% ********************* Function signature ********************************
%
% sol = vfm(csp, epsilon)
%
% *************************** Inputs **************************************
%
% csp - a well-defined "CSP" object.
%
% epsilon - the parameter for the "epsilon-net" used in SDP vector
% clustering.
%
% *************************** Outputs ************************************
%
% sol - a struct containing two fields, "folded" and "assn".
%
%   sol.folded is yet another struct, containing the fields "gurobi_output"
%   and "assn".
%
%       sol.folded.gurobi_output is the information returned from Gurobi
%       for the attempt to optimally solve the *folded* CSP.
%
%       sol.folded.assn is the assignment of variables for the *folded* CSP
%       according to the output of the Gurobi solver.
%
%   sol.assn is the assignment of variables for the input CSP after 
%   unfolding the assignment of variables for the folded CSP.
%
% *************************** BACKGROUND *********************************
%
% The Variable Folding Method is a method for generating approximate
% solutions to Constraint Satisfaction Problems that relies on the
% following 4 principles:
%
%   (1) By solving Basic SDP for "csp", we can recover a set of vectors
%   "vecs" for which the optimal SDP matrix variable is the Gram Matrix 
%   of "vecs".
%
%           For each CSP variable "i", there are |D| corresponding vectors
%           in "vecs".
%
%   (2) We can use similarily of SDP vectors (appropriately defined) to
%   lump together CSP variables. That is, if two CSP variables "i" and "j"
%   have highly similar SDP vectors, we will force CSP variables "i" and
%   "j" to take on the same value in our nearly-optimal solution.
%
%           For the Variable Folding Method, "similarity" between SDP
%           vectors is defined in the following way (to ensure provable
%           gauruntees):
%
%                   Project the SDP vectors onto a random subspace of
%                   dimension "beta". (In a proper implementation of the
%                   Variable Folding Method, beta would be extremely large.
%                   in our simplified verision of the Variable Folding
%                   Method, we simply set beta = 2.)
%
%                   Classify the projected SDP vectors over an 
%                   "epsilon-net" of the unit ball in R^beta. (In our 
%                   simplified implementation of the Variable Folding 
%                   Method, we normalize SDP vectors and and classify them
%                   over the unit sphere in R^beta).
%
%                   If two CSP variables "i" and "j" have their projected
%                   SDP vectors classified into the exact same equivalence
%                   classes, then we lump variables "i" and "j" together.
%                   If the projected SDP vectors of CSP variables "i" and
%                   "j" differ on even a single equivalence class, we do
%                   NOT merge these variables.
%
%   (3) By lumping CSP variables into a finite number of equivalence
%   classes (the maximum number of determined by beta) we define a *new*
%   CSP with a bounded number of variables. We call this new CSP the
%   "folding" of the original CSP.
%   
%   In reality, the bound on the number of variables might be extremely
%   large, but for the purposes of theoretical computer science, any
%   constant-sized bound is as good as another. Thus, the Variable Folding
%   Method states that once the new CSP is defined, one deploys an
%   EXACT solver on that CSP.
%
%   (4) The forth and final step is to take the exact solution of the
%   folded CSP (C') and "unfold" it into a solution for the original 
%   CSP (C).
%
%       If variables "i" and "j" in (C) were merged into a single variable
%       "k" in (C'), then one unfolds the solution for (C') simply by 
%       assigning variables "i" and "j" the same value as was assigned to 
%       "k" in the optimal solution to (C').
%
% ************************************************************************

%
% ************************* Implementation ********************************
%

if (nargin == 0)
    % some parameters for testing code.
    load('africaProblem_CSP_SDPL_solution_2colors.mat');
    csp = africaProblem;
    epsilon = 0.7;
end

% Step (1) : solve the SDP optimally.
sdp_sol = solve_basic_sdp(csp);

% Step (2) and (3) : classify SDP vectors & return a "folding" of the CSP. 
[~, phi_map, ~, invZ, folded_CSP] = fold(csp, sdp_sol, epsilon);
%
%   Descripton of these return variables:
%
%    for 1 <= i <= length(newvars), Z(i) is the canonical variable in the
%       original CSP associated with the i^th non-empty vector cluster.
%    for 1 <= i <= n, phi_map(i) is the canonical variable to which the i^th
%       variable in the original CSP is now mapped.
%    for 1 <= i <= n, invZ(phi_map(i)) is the label of the variable in the
%       NEW CSP to which the original variable "i" has been mapped.
%

% Step (3) (continued) : compute the optimal assignment of variables for
% the folded instance
folded_sol = solve_csp_optimally(folded_CSP);
sol.folded = folded_sol;

% Step (4) : unfolded the assignment of variables for the folded instance
% into an assignment of variables for the original instance.
n = csp.numVariables;
assignment = zeros(1,n);
for i = 1:n
    assignment(i) = folded_sol.assn(invZ(phi_map(i)));
end
sol.assn = assignment;

end