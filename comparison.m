function [Dist,Ang,H,A,G_2,V_true,m] = comparison()
% Compar the algorithm for H matrix with diffusion maps
%% Parameters:
n = 5000;
p = 0.25;
cos_alpha = 0.7; %size of the spherical cup of similar viewing angles.
%% Generating random rotations:
Rots = zeros(3,3,n);
for idx = 1:n
    Rots(:,:,idx) = randRotationMatrix;
end
%% Calculating probabilistic model and defussion maps:
[A_clean,V_true,m] = FindInvariantDistances_DiffusionMaps(Rots,cos_alpha);
[Ang,Dist] = FindOptimalAngles(Rots,A_clean);

A = Find_A_probabilistic(A_clean,p);

H_clean = zeros(n);
for Ridx = 1:n-1
    for Nidx = Ridx+1:n
        if Ang(Ridx,Nidx)~=0
            H_clean(Ridx,Nidx) = exp(1i*degtorad(Ang(Ridx,Nidx)));
            H_clean(Nidx,Ridx) = conj(H_clean(Ridx,Nidx));
        end
    end
end
H = Find_H_probabilistic(H_clean,p);

V_true_vec = reshape(V_true(1:400,1:400),1,[]);

G_1 = Find_G_matrix(H);
G_1_vec = reshape(G_1(1:400,1:400),1,[]);
figure 
scatter(V_true_vec,G_1_vec,0.05)
ylabel('G_{ij}')
xlabel('<v_i,v_j>')
title(['Calculated from H for p=',num2str(p)])

G_2 = Find_G_matrix_DiffusionMaps(A);
G_2_vec = reshape(G_2(1:400,1:400),1,[]);
figure 
scatter(V_true_vec,G_2_vec,0.05)
ylabel('G_{ij}')
xlabel('<v_i,v_j>')
title(['Calculated from A for p=',num2str(p)])
