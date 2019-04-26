function q=curve2quat_frenet(p)

% Input is a space curve p, represented as a 3xn matrix, with n being the number of sample points. Output is a curve
% q in the quaternions, obtained by lifting the Frenet framing of p. The
% quaternion is represented as a 4xn matrix.

[~,n]=size(p);

edges=polEdge(p); % Calculates derivative curve.

U=zeros(3,n);

for i=1:n
    U(:,i)=edges(:,i)/norm(edges(:,i)); % Normalized derivative curve.
end

[V,W]=frenet_frame(p);

norms=zeros(1,n);

for i=1:n
    norms(i)=norm(edges(:,i));
end

q=frame_path_to_quat(U,V,W,norms);