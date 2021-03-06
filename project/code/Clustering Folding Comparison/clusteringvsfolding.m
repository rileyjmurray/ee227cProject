clear
clc
%worldColoring;
%mycsp=africaProblem;
%mycsp=americasProblem;
epsilon=0.3:0.1:0.7;
%mycsp=read_3sat_dimacs('random_3sat_n300_c3000_1.dimacs');
pigeon_csp;
mycsp=ppcsp;

%get optimal solution
opt=solve_csp_optimally(mycsp);
temp=[evaluateObjective(mycsp,opt.assn),opt.assn];

%get clustering solution
clust_sol=solve_basic_sdp(mycsp);

    %column stack the big M matrix
    M=clust_sol.sigmaFactor;
    bigvectors=zeros(mycsp.numVariables,size(mycsp.domain,2)*size(M,1));
    for i=1:size(mycsp.domain,2)
        for j=1:mycsp.numVariables
            bigvectors(j, ((i-1)*size(M,1)+1):(i*size(M,1)))=M(:,(j-1)*size(mycsp.domain,2)+i)';
        end
    end
    for i=1:size(bigvectors,1) %normalize to unit vectors
        bigvectors(i,:)=bigvectors(i,:)/norm(bigvectors(i,:));
    end
    [x,esq,j]=kmeanlbg_solo(bigvectors,size(mycsp.domain,2)); %clustering
clust_assn=j';

%get vfm solution
vfm_assn=zeros(size(epsilon,2),mycsp.numVariables);
vfm_size=zeros(size(epsilon,2),1);
for i=1:size(epsilon,2)
    vfm_sol=vfm(mycsp,epsilon(i));
    vfm_assn(i,:)=vfm_sol.assn;
    vfm_size(i,1)=size(vfm_sol.folded.assn,2);
end

comparison=[opt.assn;clust_assn;vfm_assn];
objValue=zeros(size(comparison,1),1);
for i=1:size(comparison,1)
    objValue(i,1)=evaluateObjective(mycsp,comparison(i,:));
end

copytoexcel=[[mycsp.numVariables;mycsp.numVariables;vfm_size],objValue,comparison];