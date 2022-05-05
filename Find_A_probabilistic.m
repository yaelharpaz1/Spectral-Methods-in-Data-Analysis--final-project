function A = Find_A_probabilistic(A_clean,p)
%% Explanation: 
% Input: A_clean-size(n,n), is the indicate neighbors matrix, that
%         all the links are inferred correctly. 
%         p- probability for correct edge.
% Output: A- size(n,n), is the indicate neighbors matrix, that is a
%         small-world graph, such that with probability p for correct edge, 
%         and probability 1-p for a shortcut edge
%%
[n,~] = size(A_clean);
A = zeros(n,n);
for i_row = 1:n-1
    for j_col = i_row+1:n
        if A_clean(i_row,j_col)~=0
            x = rand;
            if x<p
                A(i_row,j_col) = A_clean(i_row,j_col);
                A(j_col,i_row) = A_clean(j_col,i_row);
            else
                idx = randi([1 n]);
                while A_clean(i_row,idx)~=0
                    idx = randi([1 n]);
                end 
                A(i_row,idx) = 1;
                A(idx,i_row) = 1;
            end
        end
    end
end

                
            