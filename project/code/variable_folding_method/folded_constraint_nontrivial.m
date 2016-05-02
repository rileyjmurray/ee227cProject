function nontriv = folded_constraint_nontrivial( constr, domain )

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

