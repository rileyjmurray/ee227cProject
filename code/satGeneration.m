function [ problem ] = satGeneration( k, m, seed )
% Generates satisfiable sat instance, i.e. OPT = 1. 
% k = num variables
% m = num constraints
% seed = (scalar) random seed for random number generator.

% set seed for repetition of experiments 
rng( seed );

% SAT relation: Scope of constraint are variables in x and y are the
% negated elements in scope. First part check if any negations are false.
% Second part checks if any positive elements evaluate true.
f = @(x,y)(sum(1-x(y))+sum(x(setdiff(1:length(x),y)))) > 0;

% generate feasible assignment which sets independently values to true with
% probability drawn from unif(0,1)
probOfTrue = rand;
feasibleAssignment = randsample(2,k,true, [ 1 - probOfTrue, probOfTrue ] ) - 1;

% Now we wish to construct a 3-SAT instance for which the above is a
% feasible solution.
constraints = cell( 1, m );
i = 1;
% Generation of random 3-SAT constraints with accept/reject generation.
while i <= m
    % generate random constraint
    % first generate how many variables are negated
    numNegated = randi( 4, 1 ) - 1;
    % select negated
    y = randsample( 3, numNegated );
    scope = randsample( k, 3, true );
    relation = @( x ) f( x , y );
    
    % construct Constraint object
    constraint = Constraint( scope, relation );
    
    % check if assignment satisfies constraint. Valid constraint if yes.
    if constraint.evaluate( feasibleAssignment )
        constraints{ i } = constraint;
        i = i + 1;
    end
end

problem = CSP( ones( 1, m ), constraints );