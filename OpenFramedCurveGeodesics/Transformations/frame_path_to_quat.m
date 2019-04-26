
% Input is a path of orthonormal frames (U,V,W) and a path in the positive
% reals r. Output is a path in quaternions, expressed as a 4xn matrix.

function q=frame_path_to_quat(U,V,W,r)

[~,n]=size(U);

rotm=zeros(3,3,n);

% Define a rotation matrix for each frame triple.
for i=1:n
    rotm(:,:,i)=[U(:,i),V(:,i),W(:,i)];
end

% Lift each rotation matrix to a unit quaternion.


for i=1:n
    q1(:,i)=rot_mat_to_quat(rotm(:,:,i));
end

% Since the unit quaternions double cover SO(3), we need to make a choice
% of lift at each point. To make the quaternionic curve continuous, we
% choose each lift to be as close as possible to the previous one. The
% fiber over a frame is +/- a quaternion.

unitq_pos=q1;
unitq_neg=-q1;

unitq(:,1)=unitq_pos(:,1);

for i=1:n-1
    if norm(unitq_pos(:,i+1)-unitq(:,i))<=norm(unitq_neg(:,i+1)-unitq(:,i))
        unitq(:,i+1)=unitq_pos(:,i+1);
    else
        unitq(:,i+1)=unitq_neg(:,i+1);
    end
end

% Define the quaternionic path by scaling each unit quaternion
% appropriately.

q=zeros(4,n);

for i=1:n
    q(:,i)=sqrt(r(i))*unitq(:,i);
end