function [X1sub, X2sub, bsub]=remove_ld_rows(X1, X2, b, tol)
%Extract a linearly independent set of rows of block matrix X = [X1, X2, b]
%
%    [X1sub, X2sub, bsub] = remove_ld_rows(X1, X2, X3, tol)

colsx1 = size(X1,2);
% colsx2 = size(X2,2); don't actually need this. keep here for readability.
X = [X1 X2 b]';
if ~nnz(X) %X has no non-zeros and hence no independent columns
    X1sub=[]; X2sub = []; bsub = [];
    return
end
if nargin < 4, tol=1e-10; end
[~, R, E] = qr(X,0); 
if ~isvector(R)
    diagr = abs(diag(R));
else
    diagr = R(1);   
end
%Rank estimation
r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
idx=sort(E(1:r));
Xsub=X(:,idx)';
X1sub = Xsub(:,1:colsx1);
X2sub = Xsub(:,(colsx1+1):(end-1));
bsub=  Xsub(:,end);
num_removed = size(X1, 1) - size(X1sub, 1);
if (num_removed > 0 )
    warning(strcat('removed _',num2str(num_removed),'_ linearly dependent constraints'));
end
       
end
       
       