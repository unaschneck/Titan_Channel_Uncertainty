function [sub_row,sub_col] = index_for_subset(row,col,start_row,start_col)
% convert row and col index into indices for row and col in a subset of the data

sub_row = row - start_row + 1;
sub_col = col - start_col + 1;


end