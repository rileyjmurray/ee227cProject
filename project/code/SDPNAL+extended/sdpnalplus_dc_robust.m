function [obj,X,s,y,S,Z,y2,v,info,runhist] =...
    sdpnalplus_dc_robust(blk,At,C,b,L,U,Bt,l,u,OPTIONS)

addpath(genpath('SDPNAL+extended/SDPT3-4.0'));
[blk,At,C,b] = sdpt3_dummy(blk,At,C,b);
rmpath(genpath('SDPNAL+extended/SDPT3-4.0'));
addpath(genpath('SDPNAL+extended/SDPNAL+v0.3'));

[obj,X,s,y,S,Z,y2,v,info,runhist] = ...
    sdpnalplus(blk,At,C,b,L,U,Bt,l,u,OPTIONS);

rmpath(genpath('SDPNAL+extended/SDPNAL+v0.3'));

end

