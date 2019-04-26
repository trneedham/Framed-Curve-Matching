function [seeds]=best_seeds_framing(p1,p2,V2,k)

% Input: A pair of space curves p1 and p2, a framing V2 of p2, and a number k.
% Output: seeds = A list of the k best shifts of p2 for minimizing geodesic 
%         distance. The geodesic distance in this step is optimized over 
%         rotations, but not reparameterizations.


[~,n]=size(p1);

% Take all shifts of p2 and compute Grassmannian distance between p1
% and each shift of p2. p1 and p2 are given Frenet framings by default. 
% We will keep the k best shifts.

Dist=zeros(1,n);

% Lift p1 to the Stiefel manifold.

q1=curve2quat_frenet(p1);
[z1,w1]=quat_to_stiefel(q1);

for j=0:n-1
    
p2j=ShiftF(p2,j);
V2j=ShiftF(V2,j);
q2j=curve2quat_arbFrame(p2j,V2j);
[z2j,w2j]=quat_to_stiefel(q2j);

distj=Grassmannian_Geodesic_Dist(z1,w1,z2j,w2j);

Dist(j+1)=distj;
end

[~,indices]=sort(Dist);

seeds=zeros(1,k);

for j=1:k
 seeds(j)=indices(j)-1;
end

