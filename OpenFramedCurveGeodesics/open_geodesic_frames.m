% Input is a pair of quaternionic curves on the L^2-sphere. Output is the 
% geodesic between the corresponding framed curves. The input curves should
% have the same number of samples and be L^2 normalized.

function [geod,frames_geod,pushoff_geod,dist]=open_geodesic_frames(q1,q2)

[~,n]=size(q1);
dt=1/n;

% Number of samples along the geodesic. This can be changed as desired.

k=100;
u=linspace(0,1,k);

% Calculate geodesic distance.

InnProd=zeros(1,n);

for i=1:n
  InnProd(i)=dot(q1(:,i),q2(:,i));
end

theta1=sum(InnProd)*dt;
theta=acos(theta1);

if (norm(theta)>.000001)
    dist=theta;
  else
    dist=0;
end

% Define the quaternionic geodesic

quat_geod=zeros(4,n,k);

for j=1:k
    quat_geod(:,:,j)=(1/sin(theta))*(sin((1-u(j))*theta)*q1(:,:)+sin(u(j)*theta)*q2(:,:));
end

% Define the geodesic between the curves and their frames

geod=zeros(3,n,k);
pushoff_geod=zeros(3,n,k);
frames_geod=zeros(3,n,k);

for j=1:k
    [geod(:,:,j),frames_geod(:,:,j),pushoff_geod(:,:,j)]=quat_to_framed_curve(quat_geod(:,:,j));
end
