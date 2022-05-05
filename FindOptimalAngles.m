function [Ang,Dist] = FindOptimalAngles(Rots,N_matrix)
%% Explanation: 
% Input: Rots-size(3,3,n), n rotation matrices in SO(3)(represents
%        projection images).
%        N_matrix- indicates neighbors. if i and j are neighbors, then
%        N_matrix(i,j)=N_matrix(j,i)=1, otherwise is 0.
% Output: Dist-size(n,n), a sparse matrix of the rotationally invariant 
%        distances, with N nonzero terms in each row. Row i represent
%        the N nearest neighbors of projection i. 
%        Ang-size(n,n), optimal alignment angles according to Dist.
%% 
[~,n] = size(N_matrix);
Dist = zeros(n,n);
Ang = zeros(n,n);
%% Find optimal angles:
fprintf('Start calculating optimal angles\n');
for Ridx = 1:n-1
    for Nidx = Ridx+1:n 
        if N_matrix(Ridx,Nidx)==1 
            R_r = ((Rots(:,:,Ridx))^(-1))*Rots(:,:,Nidx);
            C = (R_r(1,1)+R_r(2,2))/sqrt((R_r(1,1)+R_r(2,2))^2+(R_r(2,1)-R_r(1,2))^2);
            S = (R_r(2,1)-R_r(1,2))/sqrt((R_r(1,1)+R_r(2,2))^2+(R_r(2,1)-R_r(1,2))^2);
            if C >= 0
               if S >= 0   
                   angle = radtodeg(acos(C));
               else
                   angle = -radtodeg(acos(C));
               end
            else
                if S >= 0   
                    angle = 180 - radtodeg(acos(-C));
                else
                    angle = 180 + radtodeg(acos(-C));
                end
            end
            Ang(Ridx, Nidx) = angle;
            Ang(Nidx, Ridx) = mod(-angle,360);
            row = [C -S 0; S C 0; 0 0 1];
            dist = norm(Rots(:,:,Ridx)*row-Rots(:,:,Nidx), 'fro');
            Dist(Ridx, Nidx) = dist;
            Dist(Nidx, Ridx) = dist;
        end
    end
end
fprintf('Calculating optimal angles is done\n');