function [U,V,W] = rot_min_frame(p)

% For input curve p, produces a Rotation Minimizing Framing or Bishop frame
% for p. The initial frame is chosen randomly.

% p     [3,n] vector representing curve

% U     Unit tangent curve
% V     Rotation minimizing frame
% W     Completion to orthonormal basis.

[~,n] = size(p);

% Discrete tangent curve
edges=polEdge(p);

% Discrete unit tangents
for i=1:n
    U(:,i)=edges(:,i)/norm(edges(:,i));
end

V = zeros(3,n);
W = zeros(3,n);

% Initialize an initial normal vector by choosing one at random
rand_vect = rand(3,1);
tang_vect = p(:,2)-p(:,1);
proj_vect = dot(rand_vect,tang_vect)/norm(tang_vect)^2*tang_vect;
V(:,1) = (rand_vect-proj_vect)/norm(rand_vect-proj_vect);
W(:,1) = cross(U(:,1),V(:,1));

% Compute the rotation minimizing frame using the reflection method of Wang
% et. al. [Computation of Rotation Minimizing Frames, Double Reflection
% Algorithm]

for j=1:n-1
    v1 = p(:,j+1)-p(:,j);
    c1 = dot(v1,v1);
    VjL = V(:,j) - (2/c1)*dot(v1,V(:,j))*v1;
    UjL = U(:,j) - (2/c1)*dot(v1,U(:,j))*v1;
    v2 = U(:,j+1) - UjL;
    c2 = dot(v2,v2);
    V(:,j+1) = VjL - (2/c2)*dot(v2,VjL)*v2;
    W(:,j+1) = cross(U(:,j+1),V(:,j+1));
end

% Fix the RMF so that it closes
theta_diff = acos(dot(V(:,1),V(:,n)));
theta = linspace(0,theta_diff,n);
V = framing_action(p,V,theta);
W = framing_action(p,W,theta);
