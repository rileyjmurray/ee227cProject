function nontriv = folded_constraint_nontrivial( constr, domain )
%
% ********************* Function signature ********************************
%
% nontriv = folded_constraint_nontrivial( constr, domain )
%
% *************************** Inputs **************************************
%
% constr - a well-defined "Constraint" object for some CSP.
%
% domain - the domain on which all variables in the CSP (and all variables
% in this constraint) are defined.
%
% *************************** Outputs ************************************
%
% nontriv - a logical. "true" iff the constraint "constr" is non-trivial. 
% 
% A non-trivial constraint is any constraint which returns two or more 
% distinct values when looking over all possible assignments of variables
% within its scope.
%
% *************************** BACKGROUND *********************************
%
% The Variable Folding Method groups CSP variables into clusters (e.g. "i"
% and "j" into a set {i,j}) where all variables in the same cluster are
% required to take on the same value in the solution returned by the
% Variable Folding Method. 
%
% For some constraints, it is possible that all variables in its scope
% are in the same cluster. For example, given a constraint in a  max-cut
% instance with scope including variables i and j with i and j in the same
% cluster, the constraint would always evaluate to zero, regardless of the
% assignment of the cluster variable {i,j}.
%
% Trivial constraints can be removed w.l.o.g. from any CSP.
%
% ************************************************************************

num_cols = length(domain)^length(constr.scope);
lambda_space = sets2space({domain},length(constr.scope));
C =  zeros(num_cols, 1); 
for i = 1:num_cols
    tuple = linear2tuple(lambda_space, i);
    C(i) = constr.evaluate227c(tuple);
end
allsame = min(min(C) == max(C));
nontriv = 1 - allsame; 

end

