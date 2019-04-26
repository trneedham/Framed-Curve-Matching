function q=curve2quat_frenet(p)

% Input:    space curve p, represented as a 3xn matrix. 
% Output:   curve q in the quaternions, obtained by lifting the Frenet framing of p.

[~,n]=size(p);

% We need to record the parameterization speed of p. This is represented as
% a list of norms of its tangent vectors.

edges=polEdge(p);
r=zeros(1,n);
for i=1:n
    r(i)=norm(edges(:,i));
end

% Define the Frenet framing and send the weighted framing to a quaternion.
[U,V,W]=frenet_frame(p);

q=frame_path_to_quat(U,V,W,r);