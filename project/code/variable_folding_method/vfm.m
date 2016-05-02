function sol = vfm(csp, epsilon)
if (nargin == 0)
    load('africaProblem_CSP_SDPL_solution_2colors.mat');
    csp = africaProblem;
    epsilon = 0.7;
end
sdp_sol = solve_basic_sdp(csp);
[~, phi_map, ~, invZ, folded_CSP] = fold(csp, sdp_sol, epsilon);
% for 1 <= i <= length(newvars), Z(i) is the canonical variable in the
%  original CSP associated with the i^th non-empty vector cluster.
% for 1 <= i <= n, phi_map(i) is the canonical variable to which the i^th
% variable in the original CSP is now mapped.
% for 1 <= i <= n, invZ(phi_map(i)) is the label of the variable in the
% NEW CSP to which the original variable "i" has been mapped.
n = csp.numVariables;
assignment = zeros(1,n);
folded_sol = solve_csp_optimally(folded_CSP);
sol.folded = folded_sol;
for i = 1:n
    assignment(i) = folded_sol.assn(invZ(phi_map(i)));
end
sol.assn = csp.domain(assignment);

end