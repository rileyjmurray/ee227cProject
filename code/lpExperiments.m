% ----------------------------------------------------------------
% LP Experiments for 3-SAT with rounding 
% ----------------------------------------------------------------

% number of times to repeat each experiment
rep = 10;

% size of problems. Col 1 = num vars, col 2 = num constraints
expSettings = [ 
    10 20;
    10 50;
    20 40;
    20 100;
    50 50;
    50 100;
    50 200;
    ];

solveTime = - ones( size( expSettings, 1 ), rep );
roundTime = - ones( size( expSettings, 1 ), rep );
genTime = - ones( size( expSettings, 1 ), rep );
lpval = - ones( size( expSettings, 1 ), rep );
optval = - ones( size( expSettings, 1 ), rep );
roundval = -ones( size( expSettings, 1 ), rep );

for i = 1 : size( expSettings, 1 )
    for j = 1 : rep
        % set seed for repetition
        rng( expSettings( i, 1 ) * expSettings( i, 2 ) * j )
        % generate problem instance
        tic;
        % inputs: num vars, num constraints, seed
        problem = satGeneration( expSettings( i, 1 ), expSettings( i, 2 ) );
        genTime( i, j ) = toc;
        optval(i,j) = 1;
        
        tic;
        % solve LP relax
        [ muv, lambda, lpval(i,j) ] = lprelax( problem );
        solveTime( i, j ) = toc;
        
        tic;
        % apply k-SAT rounding scheme
        [ assignment ] = satRoundingLP( muv );
        roundTime( i, j ) = toc;
        roundval( i, j ) = problem.evaluateObjective( assignment );
    end
end

save( 'lpExp10rec.mat', 'roundval','optval','lpval','solveTime','genTime','roundTime' )
        
        