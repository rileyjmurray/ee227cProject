%vector clustering
%assume M (kD x kD)
k=49; %number of variables
D=4; %size of domain
%M=rand(k*D,k*D); %random M for testing
M=sol.sigmaFactor;
A=rand(2,D*size(M,1)); %random subspace
epsilon=0.5;

%Gram-Schmidt
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j=1:n
    v=A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*Q(:,j);
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end

%column stack the big M matrix
bigvectors=zeros(k,D*size(M,1));
for i=1:D
    for j=1:k
        bigvectors(j, ((i-1)*size(M,1)+1):(i*size(M,1)))=M(:,(j-1)*D+i)';
    end
end

bigvectors_proj=(Q*bigvectors')';
for i=1:size(bigvectors_proj,1)
    bigvectors_proj(i,:)=bigvectors_proj(i,:)/norm(bigvectors_proj(i,:)); %normalize to unit vectors
end

[x,esq,j]=kmeanlbg_solo(bigvectors,D) %clustering

%centers=rand(D,D*size(M_proj,1)); %random centers for classification
%centers=[eye(D);-eye(D)];
centers=epsilon_net_R2(epsilon)';
[G,J]=v_classify(bigvectors_proj,centers) %classification