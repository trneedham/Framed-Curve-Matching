% Input: Element of LS^1, entered as a list of angles, and complex curve
% (z,w), considered as an element of the Stiefel manifold.

% Output: New element of the Stiefel manifold, acted on by the elt of LS^1.

function [znew,wnew]=loops_S1_action(r,z,w)

[~,n]=size(z);

for j=1:n
    znew(j)=exp(1i*r(j))*z(j);
    wnew(j)=exp(1i*r(j))*w(j);
end

