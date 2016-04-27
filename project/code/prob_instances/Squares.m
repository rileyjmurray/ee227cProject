function [ Sproblem ] = Squares()
c = 2; %number of colors
n = 4; %size of grid
relation=@(x) x(1)~=x(2) || x(1)~=x(3) || x(1)~=x(4);

x=zeros(n,n);
x(1,:)=1:n;
for i=2:n
    x(i,:)=x(i-1,n)+1:x(i-1,n)+n;
end

%enumerate the L possibilities in the grid
possibleS = [x(1,1) x(1,n) x(n,1) x(n,n)];
increment=1;
while increment<=(n-2)
    i=1;
    while i<=(n-increment)
        j=1;
        while j<=(n-increment)
            possibleS = [possibleS; x(i,j) x(i,j+increment) x(i+increment,j) x(i+increment,j+increment)];
            j=j+1;
        end
        i=i+1;
    end
    increment=increment+1;
end

%calculate the L number
constraints=cell(size(possibleS,1),1);
for i=1:size(possibleS,1)
    constraints{i}=Constraint(possibleS(i,:), relation);
end
Sproblem=CSP(ones(size(possibleS,1),1), constraints, 0:(c-1), n^2);
end