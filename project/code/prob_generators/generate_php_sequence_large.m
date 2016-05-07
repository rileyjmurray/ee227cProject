%% test SDPs at scale with Pigeonhole Principle
pigeon_csps = cell(length(0:3:30),1);
for i = 0:3:30
    pigeon_csps{i/3+1} = pigeon_csp_fxn(10 + i, ceil((10+i)*0.8));
end
save('pigeon_csp_sequence_large.mat');