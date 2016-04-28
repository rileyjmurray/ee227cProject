function is = rowcol2svecidx(row, col)
% row = the row index at which a read / write on the matrix is desired.
% col = the column index at which a read / write on the matrix is desired.
%
% is = the index for the *vectorized* form of the matrix, where the read /
% write is desired.
temp1 = min(row,col);
temp2 = max(row,col);
row = temp1;
col = temp2;
is = (col)*(col-1)/2 + row;

end

