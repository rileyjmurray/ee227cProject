function svecx = custom_svec(x)
    s = size(x,1);
    s_bar = s * (s + 1)/2;
    svecx = zeros(1,s_bar);
    idx = 1;
    for i = 1:s
        component = x(i,i:s);
        svecx(idx:(idx + length(component) - 1)) = component;
        idx = idx + length(component);
    end
    svecx(2:(end-1)) = sqrt(2) * svecx(2:(end-1));
    svecx = sparse(svecx);
end

