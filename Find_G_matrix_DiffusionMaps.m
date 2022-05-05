function G = Find_G_matrix_DiffusionMaps(A)
%% Explanation: 
% Input: A- size(n,n), is the class averaging hermitian matrix.
% Output: G- size(n,n), is the measure of affinity between images.
%% Computing G:
fprintf('Start computing G matrix\n');
[~,n] = size(A);
D_inv = diag(1./(sum(abs(A),2)));
A_tild = D_inv*A;
[V,~] = eigs(A_tild,4);
phi = (V(:,2:4).');
phi_norm = zeros(1,n);
for i = 1:n
    phi_norm(1,i) = norm(phi(:,i),'fro');
end

G = diag(ones(1,n));
for i = 1:n-1
    for j = i+1:n 
        G(i,j) = ((phi(:,i)).'*(conj(phi(:,j))))/(phi_norm(1,i)*phi_norm(1,j));
        G(j,i) = conj(G(i,j));
    end
end
fprintf('Computing G matrix is done\n');
