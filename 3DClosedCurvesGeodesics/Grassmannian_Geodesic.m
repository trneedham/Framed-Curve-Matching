function [dist,z1new,w1new,z2new,w2new,z_geod,w_geod]=Grassmannian_Geodesic(z1,w1,z2,w2)

% Input: A pair of orthonormal bases for complex 2-planes (z1,w1),(z2,w2).
% Output: Geodesic distance between the 2-planes, 
%         Optimized bases z1new, w1new, z2new, w2new,
%         The geodesic joining the bases, z_geod, w_geod

% Define the projection matrix between the planes represented by the
% Stiefel manifold points.
A = [l2_inn_prod(z1,z2) l2_inn_prod(w1,z2);l2_inn_prod(z1,w2) l2_inn_prod(w1,w2)];

% Compute the singular value decomposition for A. This will be used to find
% optimal bases for the 2-planes and to compute geodesic distance in the
% Grassmannian.
[U,S,V] = svd(A);

z1new=V(1,1)*z1+V(2,1)*w1;
w1new=V(1,2)*z1+V(2,2)*w1;
z2new=U(1,1)*z2+U(2,1)*w2;
w2new=U(1,2)*z2+U(2,2)*w2;

z1=z1new;
w1=w1new;
z2=z2new;
w2=w2new;

Jord_z=acos(S(1,1));
Jord_w=acos(S(2,2));

% We now compute geodesic distance.
dist=sqrt(Jord_z^2+Jord_w^2);

% Normalize distance. The diameter of the Grassmannian is Sqrt(2).
if (dist>0.000001)
    dist = dist;
  else
    dist=0;
end

% Define the geodesic path parameter.
k=100;
u=linspace(0,1,k);

% Define the geodesic
for j=1:k
    z_geod(:,j)=(1/sin(Jord_z))*(sin((1-u(j))*Jord_z)*z1(:)+sin(u(j)*Jord_z)*z2(:));
    w_geod(:,j)=(1/sin(Jord_w))*(sin((1-u(j))*Jord_w)*w1(:)+sin(u(j)*Jord_w)*w2(:));
end