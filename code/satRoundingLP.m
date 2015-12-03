function [ X ] = satRoundingLP( muv )
%% LP Rounding
%
% Given the optimal muv and lambda it is easy to generate a random
% full assignment for the variables. 
%
% We have to draw based on (simply) muv; and an easy way to do so is draw
% independently for different variables. One way to do this is draw once
% from a uniform distribution and see which bin of Cumulative Distribution
% does it lie in, and assign the variable index of the bin. 
%
muv = diag(1./sum(muv,2))*muv; % Normalizing for the sake.

V = size(muv, 1);  % number of rows = number of variables

D = size(muv, 2); % Domain size for each variable

P0 = zeros(V, 1); 

P = [P0 cumsum(muv, 2)]; 
% Now each row "v" in P is the modified cdf of variable "v" 
% P(v, i) = Prob(v < i)
% We isnert the zero column, so that if the randomv ariable is less than
% muv(1) then we assign 0; and the way histcounts works given a probability
% distribution is as follows : Given P, it will count number of points in
% the given random data which lie in between the consecutive entries in P.

Z = rand(1, V); % Generate V random variables 

X = zeros(1, V); % The Full Assignment

idx = (1:V);

for v=1:V
    [Xv_idx, P(v, :)] = histcounts(Z(v), P(v, :)); 
    % Xv_idx will have exactly one non-zero entry, giving us the bin index.
                                                    
    X(v) = idx(logical(Xv_idx))-1;
    % And we assign the corresponding value (-1) to the variable.
end

