function [Dist,Ang,V_true,m] = FindInvariantDistances_rot_cup(Rots,cos_alpha)
%% Explanation: 
% Input: Rots-size(3,3,n), n rotation matrices in SO(3)(represents
%        projection images).
%        cos_alpha- a constant, where alpha is the size of the cup of similar 
%        viewing angles for every rotation matrix, chosen by the user. it 
%        implies that similar viewing angles differ by as much as alpha. 
% Output: Dist-size(n,n), a sparse matrix of the rotationally invariant 
%        distances, with N nonzero terms in each row. Row i represent
%        the N nearest neighbors of projection i. 
%        Ang-size(n,n), optimal alignment angles according to Dist.
%        V_true- size(n,n), the dot product of the true simulated viewing angles.
%        m- the number of graph edges.
%% 
%L = 360; %resolution.
[~,~,n] = size(Rots);
Dist = zeros(n,n);
Ang = zeros(n,n);
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
    %fprintf('Computing neighbors for projection %d/%d is done\n', Ridx, n);
end
fprintf('Computing neighbors for all projections is done\n');
m = nnz(N_matrix)/2;
fprintf('The number of graph edges is %d\n', m);
%%
%[x,y,z] = sphere;
%figure
%surf(x,y,z,'Facecolor','b')
%hold on
%u = Rots(:,3,1);
%scatter3(u(1,1),u(2,1),u(3,1),'filled','y')
%for i = 1:n
%    if N_matrix(1,i)~=0
%        u_tag = Rots(:,3,i);
%        scatter3(u_tag(1,1),u_tag(2,1),u_tag(3,1),'filled','r')
%    end
%end
%hold off
%title(['Example of neighborhood for cos(\alpha)=',num2str(cos_alpha)])
%% All viewing angles:
%[x,y,z] = sphere;
%figure
%surf(x,y,z,'Facecolor','b')
%hold on
%for i = 1:n
%    u = Rots(:,3,i);
%    scatter3(u(1,1),u(2,1),u(3,1),'filled','y')
%end
%hold off
%title(['All viewing angles for n=',num2str(n),' random rotations'])
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
            

           
            
            