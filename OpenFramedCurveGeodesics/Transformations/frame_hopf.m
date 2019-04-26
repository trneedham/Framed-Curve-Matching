% Takes a quaternion stored as a 4-vector [q1,q2,q3,q4] to an orthonormal frame [U,V,W]

function [U,V,W]=frame_hopf(q)

U=[q(1)^2+q(2)^2-q(3)^2-q(4)^2;2*q(2)*q(3)-2*q(1)*q(4);2*q(1)*q(3)+2*q(2)*q(4)];
V=[2*q(2)*q(3)+2*q(1)*q(4);q(1)^2-q(2)^2+q(3)^2-q(4)^2;2*q(3)*q(4)-2*q(1)*q(2)];
W=[2*q(2)*q(4)-2*q(1)*q(3);2*q(1)*q(2)+2*q(3)*q(4);q(1)^2-q(2)^2-q(3)^2+q(4)^2];
