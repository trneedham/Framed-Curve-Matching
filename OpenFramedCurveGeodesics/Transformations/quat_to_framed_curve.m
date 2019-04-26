function [p,V,pushoff]=quat_to_framed_curve(q)

% Input is a quaternionic curve, represented as a 4xn matrix. Output is a space curve p with an
% orthonormal framing V and a pushoff in the V direction.

[~,n]=size(q);

% We give a polar decomposition of q into a path of unit quaternions together with a list of norms. 
qNorms=zeros(1,n);
qUnit=zeros(4,n);

for i=1:n
    qNorms(i)=norm(q(:,i));
end

for i=1:n
    qUnit(:,i)=q(:,i)/qNorms(i);
end

% Send the unit quaternion to a path of frames.
for i=1:n
    [U(:,i),V(:,i),W(:,i)]=frame_hopf(qUnit(:,i));
end

% Squares each element in the norm list.
for i=1:n
    r(i)=qNorms(i)^2;
end

[p,pushoff]=frame_path_to_curve(U,V,W,r);