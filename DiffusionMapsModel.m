function [A,G,V_true,m] = DiffusionMapsModel()
%% Parameters:
n = 3000;
p = 0.2;
cos_alpha = 0.7; %size of the spherical cup of similar viewing angles. 
%% Generating random rotations:
Rots = zeros(3,3,n);
for idx = 1:n
    Rots(:,:,idx) = randRotationMatrix;
end
%% Calculating the dicomposition of A:
[A_clean,V_true,m] = FindInvariantDistances_DiffusionMaps(Rots,cos_alpha);

R = randn(n);

A = p*A_clean + R;
A = triu(A) + tril(A',-1); %enforce hermitian.
spec_A = eigs(A,10);
%spec_A = eig(A);
spec = sort(spec_A,'descend');

%hist_bins = 70;
%figure
%histogram(spec,hist_bins);
%title(['Histogram of the eigenvalues of A for p=',num2str(p)])
%xlabel('\lambda')

figure
bar((linspace(1,10,10)).',spec(1:10,1),1);
title(['The top 10 eigenvalues of A for p=',num2str(p)])
ylabel('\lambda')

G = Find_G_matrix_DiffusionMaps(A);
V_true_vec = reshape(V_true(1:400,1:400),1,[]);
G_vec = reshape(G(1:400,1:400),1,[]);
figure 
scatter(V_true_vec,G_vec,0.05)
ylabel('G_{ij}')
xlabel('<v_i,v_j>')
title(['p=',num2str(p)])

