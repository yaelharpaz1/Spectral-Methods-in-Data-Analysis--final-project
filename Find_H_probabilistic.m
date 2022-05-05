function H = Find_H_probabilistic(H_clean,p)
%% Explanation: 
% Input: H_clean-size(n,n), is the class averaging hermitian matrix, that
%         all the links and angles are inferred correctly. 
%         p- probability for correct edge.
% Output: H- size(n,n), is the class averaging hermitian matrix, that is a
%         small-world graph, such that with probability p for correct edge, 
%         and probability 1-p for a shortcut edge
%%
[n,~] = size(H_clean);
H = zeros(n,n);
for i_row = 1:n-1
    for j_col = i_row+1:n
        if H_clean(i_row,j_col)~=0
            x = rand;
            if x<p
                H(i_row,j_col) = H_clean(i_row,j_col);
                H(j_col,i_row) = H_clean(j_col,i_row);
            else
                idx = randi([1 n]);
                while H_clean(i_row,idx)~=0
                    idx = randi([1 n]);
                end 
                angle = randi([0 359]);
                H(i_row,idx) = exp(1i*angle);
                H(idx,i_row) = exp(1i*(-angle));
            end
        end
    end
end

                
            