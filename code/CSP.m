classdef CSP < handle
    % CSP for binary domains. It can handle non-binary domains except that
    % it does not explicitly check / able to return domains.
    
    properties
        % inputs (required and optional)
        weights
        constraints
        domain % a vector [0,...,q-1]; assume [0,1] unless otherwise stated.
        
        % calculated fields
        numConstraints 
        numVariables
        arity
    end
    
    methods
        
        function obj = CSP( weights, constraints, domain, numVariables )
            obj.domain = domain;
            obj.numVariables = numVariables;
            obj.numConstraints = length( constraints );
            % ----- Check if weights are valid -----
            if obj.numConstraints ~= length( weights )
                error('Invalid CSP: number of weights should equal the number of constraints')
            end
            if weights <= 0
                sum(weights)
                error('Invalid CSP: weights are not strictly positive')
            end
            if sum(weights) ~= 1
                disp(['normalizing weight']);
                weights= weights/sum(weights); 
            end
            obj.weights = weights;
            
            % Check if all elements in cell constraints are valid constraints
            if sum(  cellfun( @(x) isa( x, 'Constraint' ), constraints ) ) ~= obj.numConstraints
                error('Invalid CSP: constraints should contain valid Constraints objects')
            end
            obj.constraints = constraints;
            
            % determine airity
            obj.arity = 0;
            vars = [];
            for i = 1:obj.numConstraints
               tempConstr = obj.constraints{i};
               vars = union(vars,tempConstr.scope);
               obj.arity = max(obj.arity, tempConstr.arity);
            end
            obj.numVariables = length(vars);
            
        end
        
        function objectiveValue = evaluateObjective( obj, input )
            objectiveValue = 0;
            for i = 1 : obj.numConstraints
                objectiveValue = objectiveValue + obj.weights( i ) * obj.constraints{ i }.evaluate( input );
            end
        end
    end
end