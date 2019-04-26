function [p,pushoff]=frame_path_to_curve(U,V,W,r)

% Input is a path of orthonormal frames (U,V,W) and a path of positive
% reals r. Output is a path p with p'=rU and a pushoff curve 'pushoff'
% obtained as p+epsilon*V. This is a curve which is pushed off of p along the
% framing V.

% Define epsilon. This is determined by visual clarity.
epsilon=0.01;

[~,n]=size(U);
dt=1/n;

% We form the derivative of the curve by treating U as its unit tangent
% curve and r as the magnitude of its derivative.

deriv=zeros(3,n);

for i=1:n
    deriv(:,i)=r(i)*U(:,i);
end

% "Integrate" the derivative curve to obtain the curve, up to a shift in
% indices.

p1=zeros(3,n);

for i=1:3
    p1(i,:)=cumsum(deriv(i,:))*dt;
end

% Define the curve and fix the shift in indices by prepending the origin to the list.

p=zeros(3,n);

for i=1:3
    p(i,1)=0;
end

for i=2:n
    p(:,i)=p1(:,i-1);
end

pushoff=p+epsilon*V;