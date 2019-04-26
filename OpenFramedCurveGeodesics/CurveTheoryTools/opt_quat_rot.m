function q2new=opt_quat_rot(q1,q2)

% Input is a pair of quaternionic curves, given as 4xn matrices. 
% The output is a new version of the second quaternionic curve. 
% The new curve has, as a curve in C^2, been rotated by a Unitary matrix.
% The unitary matrix is chosen to minimize geodesic distance between the 
% corresponding framed space curves.

[~,n]=size(q1);

% We convert quaternionic curves to their complex representations.
z1=zeros(1,n);
w1=zeros(1,n);
z2=zeros(1,n);
w2=zeros(1,n);

for j=1:n
    z1(j)=q1(1,j)+1i*q1(2,j);
    w1(j)=q1(3,j)+1i*q1(4,j);
end

for j=1:n
    z2(j)=q2(1,j)+1i*q2(2,j);
    w2(j)=q2(3,j)+1i*q2(4,j);
end

Z1=[z1;w1];
Z2=[z2;w2];

% Alignment is performed by singular value decomposition.
A = Z1*Z2';
[U,~,V] = svd(A);
if det(A) > 0
   S = eye(2);
else
   S = eye(2);
   S(:,end) = -S(:,end);
end

Rot=U*S*V';

Z2new=zeros(2,n);

for i=1:n
    Z2new(:,i)=Rot*Z2(:,i);
end

% Convert back to quaternionic coordinates

q2new=zeros(4,n);

for j=1:n
    q2new(1,j)=real(Z2new(1,j));
    q2new(2,j)=imag(Z2new(1,j));
    q2new(3,j)=real(Z2new(2,j));
    q2new(4,j)=imag(Z2new(2,j));
end
