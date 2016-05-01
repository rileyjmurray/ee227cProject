function csp = read_sat_dimacs( filename )
    
    fileID = fopen(filename,'r');
    fgetl(fileID); % toss the first line;
    tline = fgetl(fileID); % something like "p cnf 300 3000"
    C = strsplit(tline);
    numVars = str2num(C{3});
    numClauses = str2num(C{4});
    constr = cell(numClauses,1);
    w = 1/numClauses*ones(numClauses,1);
    for i = 1:numClauses
        C = strsplit(fgetl(fileID)); % something like "-108 -213 43 0"
        v = cellfun(@str2double,C);
        relation = @(x) (x(1)*(v(1) > 0) + (1-x(1))*(v(1) < 0) || ...
            x(2)*(v(2) > 0) + (1-x(2))*(v(2) < 0) || ...
            x(3)*(v(3) > 0) + (1-x(3))*(v(3) < 0));
        constr{i} = Constraint(abs(v(1:3)), relation);
    end
    csp = CSP(w,constr,[0,1],numVars);
    fclose(fileID);
end

