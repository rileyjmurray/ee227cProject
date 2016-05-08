%% manual constraint checker : do we have all the constraints we're supposed to have?
% suppose you have a matrix variable of size 494 x 494. Well, that means
% you have 122265 scalars in that matrix variable. From the formulation of
% Basic SDP, you would expect >= 122265 linear constraints for the SDP
% representation of this SDP. 
%% load data saved at the end of solve_basic_sdp
load('/Users/RJMurray/Documents/Classes/EE227C/project/code/prob_instances/php/MACI6420160503T090319.mat');
%% grab the linear operator on the matrix variable
linOpMatVar = At{2}'; % undo the transpose done in "A>>>t<<<".
%% find the locations in the linear operator where an entry in the matrix variable is used.
for i = 1:size(linOpMatVar,1) 
    locs{i} = find(linOpMatVar(i,:)); 
end
%% assemble all of these entries into a single list, so we can see what values were used and which were not.
occured = [];
for i = 1:size(linOpMatVar,1) 
    occured = union(occured,locs{i}); % these were used 
end
didntoccur = setdiff(1:size(linOpMatVar,2),occured); % these were not
%% those indices are for the matrix variable in "svec" format. We need the indicies in the normal symmetric matrix format.
count = 1; 
for i = 1:494
    for j = i:494
        index_correspondence(count,1) = i; 
        index_correspondence(count,2) = j; 
        index_correspondence(count,3) = rowcol2svecidx(i,j); 
        count = count + 1; 
    end
end
lastIdx = count - 1;
[~, tmidx] = sort(index_correspondence(:,3));
sorted_index_correspondence = index_correspondence(tmidx,:);
%% row-col pairs of SDP matrix variable that are unconstrained.
%
% once we have this list, we need only ask *why* these are unconstrained.
% Of course, that's a little hard since "row i col j" in the SDP matrix
% variable isn't human-readable. We need to convert "row i col j" to "row
% (v,ell) col (v',ell')".
%

% the missing row and column locations
missing_row_col_locs = sorted_index_correspondence(didntoccur,1:2);
% data structures to recover the corresponding (v,ell) , (v', ell') locations
prob_vars = 1:(csp.numVariables);
dom = csp.domain;
sig_space = sets2space({prob_vars, dom});

%% inspect the locations of the missing values.
% suppose you were missing SDP matrix variable entries in row i 
% (for i \in \{3,4,...,12\}, and column 27.
%
% Why would that be? Well, run "linear2tuple(sig_space,i) and 
% linear2tuple(sig_space, 27) to see what those entries in the SDP matrix 
% variable correspond to.
%
% in our case, we found that i \in \{3,4,...,12\} corresponds to CSP
% variables 2,3,4,5,6 each taking values of 0,1. 
%
% we also found that linear2tuple(sig_space, 27) corresponded to CSP
% variable 14.
%
% we conjectued that CSP variable 14 shares no CSP constraints with CSP 
% variables 2 through 6. Indeed, after checking the scopes of all CSP 
% constraints, we did find that CSP variable 14 shares no CSP constraints
% with CSP variables 2 through 6. Thus, in our case, we had "missing" 
% constraints only because some CSP variables didn't have any CSP 
% constraints in which they interacted.
%
% linear2tuple(sig_space,3)
% linear2tuple(sig_space,4)
% linear2tuple(sig_space,5)
% linear2tuple(sig_space,6)
% linear2tuple(sig_space,7)
% linear2tuple(sig_space,8)
% linear2tuple(sig_space,9)
% linear2tuple(sig_space,10)
% linear2tuple(sig_space,11)
% linear2tuple(sig_space,12)
% linear2tuple(sig_space,27)
