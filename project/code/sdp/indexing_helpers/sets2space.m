function space = sets2space( sets, dim )
%
%       sets - a cell array of vectors. The cartesian product of these
%       vectors defines a space of tuples.
%
%       dim (optional) - if all sets are identical, then "sets" can be a
%       1x1 cell array, and "dim" can be the intended number of copies of
%       sets{1} into the rest of "sets".
%
%       space - a struct. space.offset contains indexing information that
%       needs to be computed often (and thus is best to compute once and
%       for all, and store along with the defining sets).
%       space.constituentSets is simply the cell array "sets".
%
if (nargin == 2 && length(sets) == 1)
    for i = 2:dim
        sets{i} = sets{1};
    end
end


offset = zeros(1,length(sets));
offset(end) = 1;
for i = 1:(length(sets) - 1);
    tmp = length(sets) - i;
    offset(tmp) = offset(tmp + 1) * length(sets{tmp+1});
end

space.offset = offset;
space.constituentSets = sets;

end

