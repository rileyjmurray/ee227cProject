function idx = tuple2linear( space, tuple )
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

idx = 0;
for i = 1:length(space.offset)
   loc = find(space.constituentSets{i} == tuple(i));
   if (isempty(loc))
       error('tuple not contained in space');
   end
   idx = idx + space.offset(i) * (loc - 1); 
end
idx = idx + 1;

end

