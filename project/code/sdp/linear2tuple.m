function tuple = linear2tuple( space, idx )
%
%       space - a struct containing a cell array of vectors defining the 
%       space of which "tuple" is an element, AND a vector called "offset"
%       that is used in indexing. 
%
%       tuple - a vector in domain^(dimension).
%
%       idx - a linear index according to the canonical ordering of tuples.

dim = length(space.constituentSets);
tuple = zeros(1,dim);
zidx = idx - 1;
for i = 1:dim
   tuple(i) = space.constituentSets{i}(floor(zidx / space.offset(i)) + 1);
   zidx = mod(zidx, space.offset(i)) ;
end


end

