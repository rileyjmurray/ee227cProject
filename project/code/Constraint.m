% useful relations are:
% Max cut: f = @(x) x(1) ~= x(2), where x is a 1 x 2 assignment vector
% SAT f = @(x,y)(sum(1-x(y))+sum(x(setdiff(1:length(x),y))))>0, where x is
% assignment and y is the indices of input that are negated.
classdef Constraint < handle
    properties
        scope
        arity
        relation
    end
    
    methods
        function obj = Constraint( scope, relation )
            obj.arity = length( scope );
            obj.scope = reshape(scope,1,[]);
            obj.relation = relation;
        end
        
        function value = evaluate( obj, input )
            % evaluate assignment
            inputScope = input( obj.scope );
            value = obj.relation( inputScope );
        end
    end
end