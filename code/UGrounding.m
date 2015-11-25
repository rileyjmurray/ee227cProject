%% Rounding scheme for UG problems

V = 3;
D = 2;

r = rand();

% generate 2 times D vectors with elements that are independently N(0,1)
G = normrnd( 0, 1, V * D, 2 * D );

% normalize vectors
Unorm = zeros( V *D );
for i = 1 : V * D
    Unorm( :, i ) = U( :, i ) / norm( U( :, i ) );
end

% initialize assignment
F = - ones( V, 1 );
for v = 1 : V
    xi = zeros( D, 1 );
    for l = 1 : D
        % get set size
        val = 2 * D * norm( U( :, ( v - 1 ) * D + l ) )^2;
        if val - floor( val ) <= r
            s = floor( val );
        else
            s = ceil( val );
        end
        if s == 0
            error( 'no value assigned' )
        end
        
        maxValue = - inf;
        % take inner products with gaussian and keep max value
        for i = 1 : s
            innProd = G( :, i )' * Unorm( :, ( v - 1 ) * D + l );
            if innProd > maxValue
                maxValue = innProd;
            end
        end
        xi( l ) = maxValue;
    end
    % assign F(v) to be the argmax - 1
    [ ~, argmax ] = max( xi );
    F( v ) = argmax - 1;    
    
end

    
    
    