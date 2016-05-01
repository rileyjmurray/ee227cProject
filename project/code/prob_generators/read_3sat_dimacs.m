function csp = read_3sat_dimacs( filename )
    
    fileID = fopen(filename,'r');
    ln = fgetl(fileID);
    if (ln(1) == 'c');
        ln = fgetl(fileID); % then toss that line, get a new one.
    end
    C = strsplit(ln);
    numVars = str2num(C{3});
    numClauses = str2num(C{4});
    constr = cell(numClauses,1);
    w = 1/numClauses*ones(numClauses,1);
    for i = 1:numClauses
        C = strsplit(fgetl(fileID)); % something like "-108 -213 43 0"
        v = cellfun(@str2double,C);
        if (length(v) == 4)
            relation = @(x) (x(1)*(v(1) > 0) + (1-x(1))*(v(1) < 0) || ...
                x(2)*(v(2) > 0) + (1-x(2))*(v(2) < 0) || ...
                x(3)*(v(3) > 0) + (1-x(3))*(v(3) < 0));
        elseif (length(v) == 3)
            relation = @(x) (x(1)*(v(1) > 0) + (1-x(1))*(v(1) < 0) || ...
                x(2)*(v(2) > 0) + (1-x(2))*(v(2) < 0));
        elseif (length(v) == 2)
            relation = @(x) (x(1)*(v(1) > 0) + (1-x(1))*(v(1) < 0));
            warning('have a trivial constraint (arity is 1)');
        elseif (length(v) > 4 || isempty(v))
           error('length of scope (length(v) = %d) is out of range.',length(v));
        end
        constr{i} = Constraint(abs(v(1:(length(v)-1))), relation);
    end
    csp = CSP(w,constr,[0,1],numVars);
    fclose(fileID);
end

