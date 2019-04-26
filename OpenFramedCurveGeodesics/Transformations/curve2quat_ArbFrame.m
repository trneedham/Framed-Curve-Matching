function q=curve2quat_ArbFrame(p,V)

% Input:    space curve p, normal vector field V
% Output:   curve q in the quaternions, obtained by lifting the Frenet framing of p.

[~,n]=size(p);

% We need to record the parameterization speed of p. This is represented as
% a list of norms of its tangent vectors.

edges=polEdge(p);
r=zeros(1,n);
for i=1:n
    r(i)=norm(edges(:,i));
end

for i=1:n
    U(:,i)=edges(:,i)/norm(edges(:,i));
end

for i=1:n
    W(:,i) = cross(U(:,i),V(:,i));
end

q=frame_path_to_quat(U,V,W,r);