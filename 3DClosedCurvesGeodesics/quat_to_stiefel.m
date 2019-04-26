function [z,w]=quat_to_stiefel(q)

% Input is a quaternionic curve which can be represented as a point in the
% Stiefel manifold; this means that, treated as a path in C^2, its
% coordinate functions are L^2 orthonormal.

[~,n]=size(q);

% Take the quaternionic curve to a pair of complex curves.

z=zeros(1,n);
w=zeros(1,n);

for j=1:n
    z(j)=q(1,j)+1i*q(2,j);
    w(j)=q(3,j)+1i*q(4,j);
end

% Reorthogonalize to remove any numerical error.

z=z/sqrt(l2_inn_prod(z,z));
w=w-l2_inn_prod(w,z)*z;
w=w/sqrt(l2_inn_prod(w,w));