classdef CSP < handle
    % CSP for binary domains. It can handle non-binary domains except that
    % it does not explicitly check / able to return domains.
    
    properties
        numConstraints
        weights
        constraints
    end
    
    methods
        function obj = CSP( weights, constraints )
            obj.numConstraints = length( constraints );
            % ----- Check if weights are valid -----
            if obj.numConstraints ~= length( weights )
                error('Invalid CSP: number of weights should equal the number of constraints')
            end
            if sum( weights ) ~= 1
                error('Invalid CSP: weights should sum to 1')
            end
            obj.weights = weights;
            
            % Check if all elements in cell constraints are valid constraints
            if sum(  cellfun( @(x) isa( x, 'Constraint' ), constraints ) ) ~= obj.numConstraints
                error('Invalid CSP: constraints should contain valid Constraints objects')
            end
            obj.constraints = constraints;
        end
        
        function objectiveValue = evaluateObjective( obj, input )
            objectiveValue = 0;
            for i = 1 : obj.numConstraints
                objectiveValue = objectiveValue + obj.weights( i ) * obj.constraints{ i }.evaluate( input );
            end
        end
    end
end