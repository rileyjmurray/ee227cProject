function alltuples = inctuple2alltuples( space, inctuple )

missing = find(isnan(inctuple));
alltuples = {};

if (isempty(missing))
    % This will never be the case in a recursive call. Can only be in this
    % case in the root call.
    alltuples{1} = inctuple;
elseif (length(missing) == 1)
    % This can and will happen as a result of a recursive call. 
    % There are no recusive calls after this.
    for v = space.constituentSets{missing}
       temptuple = inctuple;
       temptuple(missing) = v;
       alltuples = [alltuples, temptuple];
    end
elseif (length(missing) > 1)
   for v = space.constituentSets{missing(1)}
       inctuple(missing(1)) = v;
       sometuples = inctuple2alltuples(space, inctuple);
       alltuples = [alltuples, sometuples];
   end
end
    


end

