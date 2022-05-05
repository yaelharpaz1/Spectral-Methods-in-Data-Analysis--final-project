function G = Find_G_matrix(H)
%% Explanation: 
% Input: H- size(n,n), is the class averaging hermitian matrix.
% Output: G- size(n,n), is the measure of affinity between images.
%% Computing G:
fprintf('Start computing G matrix\n');
[~,n] = size(H);
D_inv = diag(1./(sum(abs(H),2)));
H_tild = D_inv*H;
[V,~] = eigs(H_tild,3);
psi = V.';
psi_norm = zeros(1,n);
for i = 1:n
    psi_norm(1,i) = norm(psi(:,i),'fro');
end
%psi_norm = vecnorm(psi);
G = diag(ones(1,n));
for i = 1:n-1
    for j = i+1:n 
        G(i,j) = 2*(abs(((psi(:,i)).')*(conj(psi(:,j))))/(psi_norm(1,i)*psi_norm(1,j))) - 1;
        G(j,i) = G(i,j);
    end
end
fprintf('Computing G matrix is done\n');