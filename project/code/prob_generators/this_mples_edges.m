function es = this_mples_edges( edge_space, mple )

    mple_edges = nchoosek(mple, 2);
    es = zeros(length(mple)*(length(mple)-1)/2,1);
    ct = 1;
    for i = 1:size(mple_edges, 1)
        loc_eg = mple_edges(i,:);
        es(ct) = edge_tuple2linear(edge_space,loc_eg);
        ct = ct + 1;
    end

end

