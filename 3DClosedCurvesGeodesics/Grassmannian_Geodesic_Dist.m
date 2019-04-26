function dist=Grassmannian_Geodesic_Dist(z1,w1,z2,w2)

% Input: A pair of closed space curves
% Output: Optimally aligned quaternionic curves q1new, q2new, the distance
% between the curves (dist), the geodesic joining the base curves (geod)
% and the geodesic evolution of the frames (pushoff_geod). The curves are
% by default endowed with their Frenet framings.

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
    dist = dist/sqrt(2);
  else
    dist=0;
end