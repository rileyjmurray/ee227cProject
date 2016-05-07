%% test SDPs at scale with Pigeonhole Principle
% load('pigeon_csp_sequence_large.mat');
% variable name "pigeon_csps" is a cell array.
for i = 6:length(pigeon_csps)
    sol = solve_basic_sdp(pigeon_csps{i});
    save(strcat('php_sequence_large_solution_',num2str(i),'.mat'));
end