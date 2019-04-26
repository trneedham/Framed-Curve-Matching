% Input: A pair of space curves p1, p2.
% Output: 

function [p1new,p2new,V2]=Optimal_Frame(p1,p2)

% Send the curves to quaternionic curves using their Frenet frames.
q1=curve2quat_frenet(p1);
q2=curve2quat_frenet(p2);

% Send the quaternionic curves to points in the Stiefel manifold.
[z1,w1]=quat_to_stiefel(q1);
[z2,w2]=quat_to_stiefel(q2);

% Find optimal bases for the corresponding 2-planes; that is, we find an
% optimal rotation of the framed curves.
[~,z1new,w1new,z2new,w2new,~,~]=Grassmannian_Geodesic(z1,w1,z2,w2);

[~,n]=size(q1);

% Now we find the optimal q2 under the LS^1-action. We first find the
% optimal element of LS^1.
for j=1:n
    rOpt(j)=angle(z1new(j)*conj(z2new(j))+w1new(j)*conj(w2new(j)));
end

[z2new,w2new]=loops_S1_action(rOpt,z2new,w2new);

Z1=[z1new;w1new];
Z2=[z2new;w2new];

q1new=complex_curve_to_quat(Z1);
q2new=complex_curve_to_quat(Z2);

[p1new,~,~]=quat_to_framed_curve(q1new);
[p2new,V2,~]=quat_to_framed_curve(q2new);


