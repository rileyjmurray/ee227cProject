%vector clustering
%assume M (kD x kD)
k=6; %number of variables
D=3; %size of domain
%M=rand(k*D,k*D); %random M for testing
M=sol.sigmaFactor;
A=rand(1,k*D); %random subspace

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

%project M
M_proj=Q*M;

%note: this is a column-stacking operation, but it produces row
%vectors because that is what the clustering function expects
vectors=zeros(k,D*size(M_proj,1));
for i=1:D
    for j=1:k
        %vectors(j, ((i-1)*k*D+1):(i*k*D))=M(:,(j-1)*D+i)'; %unprojected
        vectors(j, ((i-1)*size(M_proj,1)+1):(i*size(M_proj,1)))=M_proj(:,(j-1)*D+i)'; %projected
    end
end
for i=1:size(vectors,1)
    vectors(i,:)=vectors(i,:)/norm(vectors(i,:)); %normalize to unit vectors
end

[x,esq,j]=kmeanlbg_solo(vectors,D) %clustering

%centers=rand(D,D*size(M_proj,1)); %random centers for classification
centers=eye(D);
[G,J]=v_classify(vectors,centers) %classification