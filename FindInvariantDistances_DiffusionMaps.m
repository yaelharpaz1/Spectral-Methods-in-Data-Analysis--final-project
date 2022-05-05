function [A,V_true,m] = FindInvariantDistances_DiffusionMaps(Rots,cos_alpha)
%% Explanation: 
% Input: Rots-size(3,3,n), n rotation matrices in SO(3)(represents
%        projection images).
%        cos_alpha- a constant, where alpha is the size of the cup of similar 
%        viewing angles for every rotation matrix, chosen by the user. it 
%        implies that similar viewing angles differ by as much as alpha. 
% Output: A- indicates neighbors. if i and j are neighbors, then
%        A(i,j)=A(j,i)=1, otherwise is 0.
%        V_true- size(n,n), the dot product of the true simulated viewing angles.
%        m- the number of graph edges.
%%
[~,~,n] = size(Rots);
V_true = diag(ones(1,n)); % dot product of the true simulated viewing angles.
N_matrix = zeros(n,n); % indicates neighbors.
%% Find neighbors of every rotation matrix: 
fprintf('Start computing neighbors\n');
for Ridx = 1:n-1
    for Nidx = Ridx+1:n
        V_true(Ridx,Nidx) = (Rots(:,3,Ridx).')*Rots(:,3,Nidx);
        V_true(Nidx,Ridx) = V_true(Ridx,Nidx);
        if V_true(Ridx,Nidx) > cos_alpha
            N_matrix(Ridx,Nidx) = 1;
            N_matrix(Nidx,Ridx) = 1;
        end
    end
end
fprintf('Computing neighbors for all projections is done\n');
m = nnz(N_matrix)/2;
fprintf('The number of graph edges is %d\n', m);
A = N_matrix;
