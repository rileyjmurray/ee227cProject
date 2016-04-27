[obj,X,~,~,~,~,~,~,info,runhist] = ...
    sdpnalplus_dc_robust(blk,At,C,b,L,U,Bt,l,u,OPTIONS);

[L,D] = ldl(X{2});
M = L*sqrt(D);

sol.lambdaVar = X{1};
sol.sigmaVar = X{2};
sol.sigmaFactor = M;
sol.obj = obj;
sol.termcode = info.termcode;
sol.all_info = info;
sol.run_hist = runhist;

save('ramsey4n18.mat');