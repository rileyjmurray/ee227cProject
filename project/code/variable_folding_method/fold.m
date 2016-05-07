function [A_map, phi_map, Z, invZ, folded_CSP] = fold(csp, sdpsol, epsilon, seed)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin == 0)
        load('africaProblem_CSP_SDPL_solution_2colors.mat');
        csp = africaProblem;
        epsilon = 0.7;
        rng(0);
    elseif (nargin == 3)
        rng(0);
    else
        rng(seed); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% random projection %%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = 2; %... don't change this.
    n = csp.numVariables;
    D = length(csp.domain);
    PHI_proj = normc(normrnd(0,1/beta,beta,n*D));
    uvecs = PHI_proj * sdpsol.sigmaFactor;
    plotv(uvecs);
    % ... technically, supposed to identify and drop the "bad" constraints
    % given this projection.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% discretization %%%%%%%%%%%%%%%%%%%%%%%%%%
    [A_map, phi_map, Z, invZ] = vfm_mappings(uvecs, n, D, epsilon);
    % for 1 <= i <= length(newvars), Z(i) is the canonical variable in the
    %  original CSP associated with the i^th non-empty vector cluster.
    % for 1 <= i <= n, phi_map(i) is the canonical variable to which the i^th
    % variable in the original CSP is now mapped.
    % for 1 <= i <= n, invZ(phi_map(i)) is the label of the variable in the
    % NEW CSP to which the original variable "i" has been mapped.
    
    %%%%%%%%%%%%%%%%%%%%%%%% construct folded CSP %%%%%%%%%%%%%%%%%%%%%%%%%
    count = 0;
    keepers = [];
    clear folded_constraints;
    for j = 1:csp.numConstraints
        C = csp.constraints{j};
        ns = zeros(1,length(C.scope));
        for i = 1:length(ns)
            oldvar = C.scope(i);
            folded_oldvar = phi_map(oldvar);
            newvar = invZ(folded_oldvar);
            % display([j, oldvar, folded_oldvar, newvar]);
            ns(i) = newvar;
        end
        [newscope,~,icount] = unique(ns,'stable');
        newrelation = @(x) C.relation(x(icount));
        newConstraint = Constraint(newscope, newrelation);
        if (folded_constraint_nontrivial(newConstraint, csp.domain))
            keepers(end+1) = j;
            folded_constraints{count+1} = newConstraint;
            count = count + 1;
        end
    end
    folded_CSP = CSP(csp.weights(keepers), folded_constraints, csp.domain, length(Z));
end
 