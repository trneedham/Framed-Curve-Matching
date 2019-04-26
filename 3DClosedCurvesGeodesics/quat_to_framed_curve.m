function [p,V,pushoff]=quat_to_framed_curve(q)

% Input is a quaternionic curve, represented as a 4xn matrix. Output is a space curve p with an
% orthonormal framing V and a pushoff in the V direction.

[U,V,W,r]=quat_to_frame_curve(q);
[p,pushoff]=frame_path_to_curve(U,V,W,r);