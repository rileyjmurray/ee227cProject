%% summary
% write to DIMACS format as n-SAT (... lol)
% convert n-SAT to 3-SAT with the python function borrowed from
% StackExchange.
% load the 3-SAT representation in as a CSP

%% params
n = 15; % number of items
m = 11; % number of boxes
numlits = n*m; % variables are blocked by box (i.e. x_{i1} is in a block, x_{i2} is in the next block)
numclauses = n + m*(n*(n-1)/2);
% normally, x_{ij} says "item i in box j".
%% write opening lines
filename = strcat('pigeon_nSAT_n',num2str(n),'_m',num2str(m),'.dimacs');
fileID = fopen(filename,'w');
fprintf(fileID, 'p cnf %d %d \n', numlits, numclauses);

%%
insomebox = cell(n,1);
for i = 1:n
    insomebox{i} = i:n:numlits; % item "i" is in some box.
    fprintf(fileID, ' %d', insomebox{i});
    fprintf(fileID, ' %d \n', 0);
end
%%
notbothinbox = cell(m*(n*(n-1)/2),1);
c = 0;
for j = 1:m
   for i1 = 1:n
       for i2 = (i1+1):n
           i1_in_box_j = (j-1)*n + i1;
           i2_in_box_j = (j-1)*n + i2;
           notbothinbox{c+1} = [-i1_in_box_j, -i2_in_box_j];
           fprintf(fileID, ' %d', notbothinbox{c+1});
           fprintf(fileID, ' %d \n', 0);
           c = c + 1;
       end
   end
end
if (c ~= length(notbothinbox))
    error('something went wrong');
end
fclose(fileID);
%% convert
commandStr = char(strcat(...
    {'python prob_generators/convert_diamcs_sat_to_dimacs_3sat.py'},...
    {' '},{filename}));
[status, commandOut] = system(commandStr);
if status~=0
    warning('something went wrong');
end
%% read in as CSP
ppcsp = read_3sat_dimacs(char(strcat(filename,'_3sat')));

