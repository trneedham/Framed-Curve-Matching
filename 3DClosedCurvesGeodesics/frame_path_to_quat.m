function q=frame_path_to_quat(U,V,W,r)

% Input is a path of orthonormal frames (U,V,W) and a path in the positive
% reals r. Output is a path in quaternions, expressed as a 4xn matrix.

[~,n]=size(U);

for i=1:n
    rotm(:,:,i)=[U(:,i),V(:,i),W(:,i)];
end

for i=1:n
    q1(i,:)=rot_mat_to_quat(rotm(:,:,i));
end

% Since the unit quaternions double cover SO(3), we need to make a choice
% of lift at each point. To make the quaternionic curve continuous, we
% choose each lift to be as close as possible to the previous one. The
% fiber over a frame is +/- a quaternion.

% Note: If the input curve has corners, then this choice of quaternionic
% path could cause problems!

unitq_pos=q1';
unitq_neg=-q1';

unitq(:,1)=unitq_pos(:,1);

for i=1:n-1
    if norm(unitq_pos(:,i+1)-unitq(:,i))<=norm(unitq_neg(:,i+1)-unitq(:,i))
        unitq(:,i+1)=unitq_pos(:,i+1);
    else
        unitq(:,i+1)=unitq_neg(:,i+1);
    end
end

for i=1:n
    q(:,i)=sqrt(r(i))*unitq(:,i);
end