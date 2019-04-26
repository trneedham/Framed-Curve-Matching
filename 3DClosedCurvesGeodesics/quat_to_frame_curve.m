function [U,V,W,r]=quat_to_frame_curve(q)

% Input is a quaternionic curve as a 4xn matrix. Output is a path of orthonormal frames
% [U,V,W] and a path of positive real numbers r.

[~,n]=size(q);

% Decompose the quaternionic path into polar coordinates.
for i=1:n
    qNorms(i)=norm(q(:,i));
end

for i=1:n
    qUnit(:,i)=q(:,i)/qNorms(i);
end

for i=1:n
    [U(:,i),V(:,i),W(:,i)]=frame_hopf(qUnit(:,i));
end

for i=1:n
    r(i)=qNorms(i)^2;
end
