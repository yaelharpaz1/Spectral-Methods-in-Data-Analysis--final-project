function [Dist,Ang,H,G,V_true,m] = ProbabilisticModel()
%% Parameters:
n = 3000;
p = 1;
cos_alpha = 0.7; %size of the spherical cup of similar viewing angles. 
%% Generating random rotations:
Rots = zeros(3,3,n);
for idx = 1:n
    Rots(:,:,idx) = randRotationMatrix;
end
%% Calculating the dicomposition of H:
[Dist,Ang,V_true,m] = FindInvariantDistances_rot_cup(Rots,cos_alpha);
fprintf('Start computing H matrix\n');
H_clean = zeros(n);
for Ridx = 1:n-1
    for Nidx = Ridx+1:n
        if Ang(Ridx,Nidx)~=0
            H_clean(Ridx,Nidx) = exp(1i*degtorad(Ang(Ridx,Nidx)));
            H_clean(Nidx,Ridx) = conj(H_clean(Ridx,Nidx));
        end
    end
end

R = randn(n);
%H = Find_H_probabilistic(H_clean,p);

H = p*H_clean + R;
H = triu(H) + tril(H',-1); %enforce hermitian.
fprintf('Computing H matrix is done\n');
spec_H = eigs(H,50);
%spec_H = eig(H);
spec = sort(spec_H,'descend');

%hist_bins = 70;
%figure
%histogram(spec,hist_bins);
%title(['Histogram of the eigenvalues of H for p=',num2str(p)])
%xlabel('\lambda')

figure
bar((linspace(1,50,50)).',spec(1:50,1),1);
title(['The top 50 eigenvalues of H for p=',num2str(p)])
ylabel('\lambda')

G = Find_G_matrix(H);
V_true_vec = reshape(V_true(1:400,1:400),1,[]);
G_vec = reshape(G(1:400,1:400),1,[]);
figure 
scatter(V_true_vec,G_vec,0.05)
ylabel('G_{ij}')
xlabel('<v_i,v_j>')
title(['p=',num2str(p)])

