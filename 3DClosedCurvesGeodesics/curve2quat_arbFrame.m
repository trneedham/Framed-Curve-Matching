function q=curve2quat_arbFrame(p,V)

% Input is a space curve p, represented as a 3xn matrix and a framing V. Output is a curve
% q in the quaternions, obtained by lifting the framed curve (p,V).

[d,n]=size(p);

edges=polEdge(p);

for j=1:n
    U(:,j)=edges(:,j)/norm(edges(:,j));
    W(:,j)=cross(U(:,j),V(:,j));
end

norms=zeros(1,n);

for i=1:n
    norms(i)=norm(edges(:,i));
end

q=frame_path_to_quat(U,V,W,norms);