function tuple = edge_linear2tuple( space, idx )
%
%       space - a struct containing a cell array of vectors defining the 
%       space of which "tuple" is an element, AND a vector called "offset"
%       that is used in indexing. 
%
%       tuple - a vector in domain^(dimension).
%
%       idx - a linear index according to the canonical ordering of tuples.

error('don''t call me! Pre-compute this once and for all with "nchoosek(1:n,2)"');

end

