function [ropt,q2opt]=optimal_framing(q1,q2)

% Input is a pair of quaternionic curves corresponding to a pair framed curves.
% The output is a new q2, optimized over the paths S1 action by frame
% twisting. The input quaternions should lie on the L^2-sphere.

[~,n]=size(q1);

% Now we find the optimal q2 under the PS^1-action. This is done by a
% formula found in the paper corresponding to this program. The process
% involves converting the quaternionic curve into its complex
% representation.
z1=zeros(1,n);
w1=zeros(1,n);
z2=zeros(1,n);
w2=zeros(1,n);

for j=1:n
    z1(j)=q1(1,j)+1i*q1(2,j);
    w1(j)=q1(3,j)+1i*q1(4,j);
end

for j=1:n
    z2(j)=q2(1,j)+1i*q2(2,j);
    w2(j)=q2(3,j)+1i*q2(4,j);
end

ropt=zeros(1,n);

for j=1:n
    ropt(j)=angle(z1(j)*conj(z2(j))+w1(j)*conj(w2(j)));
end

% Let the optimal rotation act on q2
for j=1:n
    z2(j)=exp(1i*ropt(j))*z2(j);
    w2(j)=exp(1i*ropt(j))*w2(j);
end

q2opt=zeros(4,n);

% Convert back to quaternionic representation
for j=1:n
    q2opt(1,j)=real(z2(j));
    q2opt(2,j)=imag(z2(j));
    q2opt(3,j)=real(w2(j));
    q2opt(4,j)=imag(w2(j));
end
