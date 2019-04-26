% Input is a quaternionic curve, given as a 4xn matrix---the rows are
% coordinate functions and each column is a sample along the curve. The
% output is the same curve, written as a pair of curves in the complex
% numbers.

function [z,w]=quat_curve_to_complex(q);

[d,n]=size(q);

for j=1:n
    z(j)=q(1,j)+i*q(2,j);
    w(j)=q(3,j)+i*q(4,j);
end