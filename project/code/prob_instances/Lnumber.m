function [ Lproblem ] = Lnumber()
c = 2; %number of colors
n = 4; %size of grid
x=zeros(n,n);

    function [bool] = Lfunction (x)
        %enumerate the L possibilities in the grid
        en=size(x,1);
        possibleL = [x(1,1) x(en,1) x(en,en)];
        increment=1;
        while increment<=(en-2)
            i=1;
            while i<=(en-increment)
                j=1;
                while j<=(en-increment)
                    possibleL = [possibleL; x(i,j) x(i+increment,j) x(i+increment,j+increment)];
                    j=j+1;
                end
                i=i+1;
            end
            increment=increment+1;
        end
        
        %calculate the L number
        L=0;
        for i=1:size(possibleL,1)
            L=L+(possibleL(i,1)==possibleL(i,2) && possibleL(i,2)==possibleL(i,3));
        end
        if L==0
            bool=1;
        else
            bool=0;
        end
    end

constraints=cell(1,1);
constraints{1}=Constraint(x,@Lfunction);
Lproblem=CSP(1,constraints,0:(c-1),n^2);
end