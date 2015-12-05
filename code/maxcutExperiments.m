% ----------------------------------------------------------------
% Experiments for Max-CUT with rounding 
% ----------------------------------------------------------------

% number of times to repeat each experiment
rep = 5;

% size of problems. Col 1 = num vars, col 2 = num constraints
expSettings = [ 
    5 4;
    10 50;
    20 40;
%     20 100;
%     50 50;
%     50 100;
%     50 200;
%     100 100;
%     100 200;
%     100 400;
    ];
genTime = - ones( size( expSettings, 1 ), rep );
optval = - ones( size( expSettings, 1 ), rep );

ugsolveTime = - ones( size( expSettings, 1 ), rep );
ugroundTime = - ones( size( expSettings, 1 ), rep );
ugsdpval = - ones( size( expSettings, 1 ), rep );
ugroundval = -ones( size( expSettings, 1 ), rep );

lpsolveTime = - ones( size( expSettings, 1 ), rep );
lproundTime = - ones( size( expSettings, 1 ), rep );
lpval = - ones( size( expSettings, 1 ), rep );
lproundval = -ones( size( expSettings, 1 ), rep );

gwsolveTime = - ones( size( expSettings, 1 ), rep );
gwroundTime = - ones( size( expSettings, 1 ), rep );
gwsdpval = - ones( size( expSettings, 1 ), rep );
gwroundval = -ones( size( expSettings, 1 ), rep );


for i = 1 : size( expSettings, 1 )
    for j = 1 : rep
        % set seed for repetition
        rng( expSettings( i, 1 ) * expSettings( i, 2 ) * j )
        % generate problem instance
        tic;
        % inputs: num vars, num constraints, seed
        V =  expSettings( i, 1 );
        D =  2;
        C = expSettings( i, 2 ) ;
        problem = maxcutBipartite(V, C);
        genTime( i, j ) = toc;
        optval(i,j) = 1;
        
        tic;
        
        % Solve different SDP GW
        [sigmaVargw, gwsdpval(i, j)] = gwSDP(problem);
        gwsolveTime(i, j) = toc;
        
        tic;
        % solve SDP relax UG
        [ sigmaVarug, lambdasdp, ugsdpval(i,j) ] = constructAndSolveSDP( problem );
        ugsolveTime( i, j ) = toc;
        
        tic;
        %solve lp relax
        [ muv, lambdalp, lpval(i,j) ] = lprelax( problem );
        lpsolveTime( i, j ) = toc;
        
        
        tic;
        % apply gw rounding scheme
        [assignment] = gwRound(sigmaVargw);
        gwroundTime( i, j ) = toc;
        gwroundval( i, j ) = problem.evaluateObjective( assignment );
        
        tic;
        % apply UG rounding scheme
        [ assignment ] = UGrounding(sigmaVarug, V, D);
        ugroundTime( i, j ) = toc;
        ugroundval( i, j ) = problem.evaluateObjective( assignment );
        
        tic;
        % apply LP rounding scheme
        [ assignment ] = satRoundingLP(muv);
        lproundTime( i, j ) = toc;
        lproundval( i, j ) = problem.evaluateObjective( assignment );
        
    end
end

save( ['maxcut' num2str( rep ) 'rec.mat']) ; %, 'roundval','optval','sdpval','solveTime','genTime','roundTime' )
        
        