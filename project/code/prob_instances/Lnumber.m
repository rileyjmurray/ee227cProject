function Lproblem = Lnumber()

c = 2; %number of colors
n = 4; %size of grid
relation=@(x) x(1)~=x(2) || x(1)~=x(3);

x=zeros(n,n);
x(1,:)=1:n;
for i=2:n
    x(i,:)=x(i-1,n)+1:x(i-1,n)+n;
end

%enumerate the L possibilities in the grid
possibleL = [x(1,1) x(n,1) x(n,n)];
increment=1;
while increment<=(n-2)
    i=1;
    while i<=(n-increment)
        j=1;
        while j<=(n-increment)
            possibleL = [possibleL; x(i,j) x(i+increment,j) x(i+increment,j+increment)];
            j=j+1;
        end
        i=i+1;
    end
    increment=increment+1;
end

%calculate the L number
constraints=cell(size(possibleL,1),1);
for i=1:size(possibleL,1)
    constraints{i}=Constraint(possibleL(i,:), relation);
end
Lproblem=CSP(ones(size(possibleL,1),1), constraints, 0:(c-1), n^2);
a = 1;
end