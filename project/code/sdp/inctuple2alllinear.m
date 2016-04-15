function indices = inctuple2alllinear( space, inctuple )

alltuples = inctuple2alltuples(space, inctuple);
indices = zeros(1,length(alltuples));
for i = 1:length(alltuples)
    indices(i) = tuple2linear(space, alltuples{i});
end


end

