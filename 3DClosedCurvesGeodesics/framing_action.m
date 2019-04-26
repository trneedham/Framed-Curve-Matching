function Vnew = framing_action(p,V,theta)

% Input: base curve p, initial framing V, angle function theta
% Output: New frame obtained by acting on V according to the function theta

[~,n] = size(p);

% Discrete tangent curve
edges = polEdge(p);

% Discrete unit tangents
for i=1:n
    U(:,i) = edges(:,i)/norm(edges(:,i));
end

% Complete the frame
for i=1:n
    W(:,i) = cross(U(:,i),V(:,i));
end

% Act on the frame
for i=1:n
    Vnew(:,i) = V(:,i)*cos(theta(i))+W(:,i)*sin(theta(i));
end

