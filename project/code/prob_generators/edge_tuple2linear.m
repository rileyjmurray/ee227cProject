function idx = edge_tuple2linear( space, tuple )
%
%       space - a struct containing a cell array of vectors defining the 
%       space of which "tuple" is an element, AND a vector called "offset"
%       that is used in indexing. 
%
%       tuple - [t1, t2, ... t(length(space))], with t(i) \in space{i}
%
%       idx - a linear index according to the canonical ordering of tuples.
%       (for two elements ti and tj in "tuple", if i > j then i is "less
%       strided" than j).
%
n = space.constituentSets{1}(end); % assumes "space" is 1:n
if (min(tuple) ~= tuple(1) || max(tuple) ~= tuple(2))
   warning('edges must have the smaller vertex as the first entry. ');
   warning('redefining input tuple to be consistent with edge ordering rules');
   tuple = [min(tuple), max(tuple)];
end
idx = tuple2linear(space, tuple);
idx = idx - tuple(1)*(tuple(1)+1)/2;


end

