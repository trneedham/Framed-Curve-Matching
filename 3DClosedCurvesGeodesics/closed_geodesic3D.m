function [dist,geod,vect_geod,pushoff_geod]=closed_geodesic3D(p1,p2)

% Input: A pair of closed space curves
% Output: Optimally aligned quaternionic curves q1new, q2new, the distance
% between the curves (dist), the geodesic joining the base curves (geod)
% and the geodesic evolution of the frames (pushoff_geod). The curves are
% by default endowed with their Frenet framings.

    N=200;

    p1 = ReSampleCurve(p1,N);
    p2 = ReSampleCurve(p2,N);
    
% Preprocess the curves to have the same length.
    p1=preprocess_curve(p1);
    p2=preprocess_curve(p2);
    
    
% Lift the Frenet framed curves to quaternions.
q1=curve2quat_frenet(p1);
q2=curve2quat_frenet(p2);

% Send the quaternionic curves to their Stiefel manifold representations.
[z1,w1]=quat_to_stiefel(q1);
[z2,w2]=quat_to_stiefel(q2);

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
    dist = 2*dist/(sqrt(2)*pi);
  else
    dist=0;
end

% Define the geodesic path parameter.
k=100;
u=linspace(0,1,k);

Z1new=[z1new;w1new];
Z2new=[z2new;w2new];

% Send the optimal bases back to quaternionic curves.
q1new=complex_curve_to_quat(Z1new);
q2new=complex_curve_to_quat(Z2new);

% Define the geodesic
for j=1:k
    z_geod(:,j)=(1/sin(Jord_z))*(sin((1-u(j))*Jord_z)*z1(:)+sin(u(j)*Jord_z)*z2(:));
    w_geod(:,j)=(1/sin(Jord_w))*(sin((1-u(j))*Jord_w)*w1(:)+sin(u(j)*Jord_w)*w2(:));
end

for j=1:k
    Z_geod(1,:,j)=z_geod(:,j);
    Z_geod(2,:,j)=w_geod(:,j);
    quat_geod(:,:,j)=complex_curve_to_quat(Z_geod(:,:,j));
end

for j=1:k
    [geod(:,:,j),vect_geod(:,:,j),pushoff_geod(:,:,j)]=quat_to_framed_curve(quat_geod(:,:,j));
end

for j=1:k
    geod(:,end,j)=[0;0;0];
end

for j=1:k
    vect_geod(:,end,j)=vect_geod(:,1,j);
end
