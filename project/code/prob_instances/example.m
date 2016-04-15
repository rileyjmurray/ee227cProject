% EXAMPLE
% want to create max SAT csp with variables x1, x2, x3
% Constraints: 
% 1) x1 or not x2 or x3, weight = .6
% 2) x2 or not x3, weight = .4

% set relation, x is an assignment (only for variables in scope) and y are
% the negated variables in constraint
r = @(x,y)(sum(1-x(y))+sum(x(setdiff(1:length(x),y))))>0;

% define constraint 1). @(x) r( x, [ 2] ) is an anonymous function that
% takes as input the assignment (for all variables) the values in x,
% where 1st and 3rd variable need to be true and the 2nd false.
C1 = Constraint( [ 1, 2, 3 ], @(x) r( x, [ 2 ] ) );

% evaluating [ 0 1 0] --> x1 = 0, x2 = 1, x3 = 0 correctly returns 0
C1.evaluate( [ 0 1 0 ] )
% evaluating [ 1 0 0] --> x1 = 1, x2 = 0, x3 = 0 correctly returns 1
C1.evaluate( [ 1 0 0 ] )
% accessing its scope gives [ 1 2 3], i.e. the variables are x1, x2, x3
C1.scope

% similarly we have 2) (the [2] here says that the second variable in the
% scope should be negated )
C2 = Constraint( [ 2 3], @(x) r( x, [2] ) );
% evaluating [ 0 0 1 ] for 2) --> x1= 0, x2 = 0, x3 = 1 correctly returns 0
% (note that the input to the evaluate functions is the whole assignment
% and NOT the variables in its scope.
C2.evaluate( [ 0 0 1] )

% set of constraints:
constraints = { C1, C2 };
% and weights
weights = [ .6 .4 ];

% construct CSP
problem = CSP( weights, constraints );
% to determine the objective value for a given assignment x1 = 0, x2 = 0,
% x3 = 1. Observe that we get .6 as expected:
problem.evaluateObjective( [ 0 0 1 ] )

