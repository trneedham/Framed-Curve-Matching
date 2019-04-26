% The input for this function is a space curve, given as a list of points
% parameterized by time. The output is the discrete tangent curve.
% The code addresses the issue that the discrete derivative decreases the
% number of samples by 1.

% This code assumes that the curve is parameterized on the interval [0,1].

function ed = polEdge(p)

[n,T] = size(p);

dt=1/T;

for i=1:n
  polEdge(i,:)=diff(p(i,:))/dt;
end

% If p is a closed curve then we want its derivative curve to also be
% closed. We define the failure to close number:
closure_diff=norm(p(:,1)-p(:,T));

% If the failure to close is small, then we force the curve to close. If
% not we extend the derivative curve so that our indexing doesn't
% change. The extension will be in the direction of the previous tangent
% vector (so it will have a small straight piece at the end).
if closure_diff<.001
    for i=1:n
        polEdge(i,T)=polEdge(i,1);
    end
else
    polEdge(:,T)=polEdge(:,T-1)+polEdge(:,T-1)-polEdge(:,T-2);
end

ed=polEdge;