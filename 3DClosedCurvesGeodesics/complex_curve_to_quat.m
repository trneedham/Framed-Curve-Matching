% Input is a curve in C^2, given as a 2xn complex matrix. Output
% is the same curve, expressed as a curve in R^4 (or the quaternions).

function q=complex_curve_to_quat(Z)

[~,n]=size(Z);

for j=1:n
    q(1,j)=real(Z(1,j));
    q(2,j)=imag(Z(1,j));
    q(3,j)=real(Z(2,j));
    q(4,j)=imag(Z(2,j));
end

