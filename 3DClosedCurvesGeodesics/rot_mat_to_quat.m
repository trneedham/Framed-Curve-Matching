function q = rot_mat_to_quat(R)

% This is a modification of rotm2q from the "quaternion" package of Reichlin and Carbajal.
% It was originally written for octave.
% The modification is made so that the particular version of the 
% The Hopf map in our conventions is a left inverse.
% The input for this function is a rotation matrix (i.e. an element of SO(3))
% and the output is a quaternion. 

  if (nargin ~= 1)
    print_usage ();
  end

  T = trace (R);
  
  if (T > 0)
    s = 0.5 / sqrt (T+1);
    w = 0.25 / s;
    x = (R(2,3) - R(3,2)) * s;
    y = (-R(1,3) + R(3,1)) * s;
    z = (-R(2,1) + R(1,2)) * s;
  else

    if (R(1,1) > R(2,2) && R(1,1) > R(3,3))
      s = 2 * sqrt (1 + R(1,1) - R(2,2) - R(3,3));
      w = (-R(3,2) + R(2,3)) / s;
      x = 0.25 * s;
      y = (R(1,2) + R(2,1)) / s;
      z = (R(1,3) + R(3,1)) / s;
    elseif (R(2,2) > R(3,3))
      s = 2 * sqrt (1 + R(2,2) - R(1,1) - R(3,3));
      w = (-R(1,3) + R(3,1)) / s;
      x = (R(1,2) + R(2,1)) / s;
      y = 0.25 * s;
      z = (R(2,3) + R(3,2) ) / s;
    else
      s = 2 * sqrt (1 + R(3,3) - R(1,1) - R(2,2));
      w = (-R(2,1) + R(1,2)) / s;
      x = (R(1,3) + R(3,1)) / s;
      y = (R(2,3) + R(3,2)) / s;
      z = 0.25 * s;
    end
  end

  q = [w, x, y, z];



