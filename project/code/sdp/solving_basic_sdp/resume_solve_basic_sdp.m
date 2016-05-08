% Use this script when you have loaded input data for SDPNAL+, but never
% ran SDPNAL+ itself.

[obj,X,~,~,~,~,~,~,info,runhist] = ...
    sdpnalplus(blk,At,C,b,L,U,Bt,l,u,OPTIONS);

[L,D] = ldl(X{2});
M = L*sqrt(D);

sol.lambdaVar = X{1};
sol.sigmaVar = X{2};
sol.sigmaFactor = M;
sol.obj = obj;
sol.termcode = info.termcode;
sol.all_info = info;
sol.run_hist = runhist;

save('resumed_solving_some_Basic_SDP.mat');