function edges = polEdge(p)

% The input for this function is a (open or closed) space curve, given as a 
% 3xn matrix with n samples denoting parameter value. The output is the discrete tangent curve.

[n,T] = size(p);

dt=1/T;

for i=1:n
  edges(i,:)=diff(p(i,:))/dt;
end

% If this curve was closed, this will fix the tangent indicatrix to be
% closed as well.
if p(:,1)-p(:,T)==0 
    edges(:,T)=edges(:,1);
else
    edges(:,T)=edges(:,T-1);
end
