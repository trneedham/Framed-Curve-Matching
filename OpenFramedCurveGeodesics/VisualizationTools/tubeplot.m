function [x,y,z]=tubeplot(curve,frame,r,n)

% This code is adapted from that of Janus H. Wesenberg, july 2004. The
% original code used Frenet framing. This has been updated to use an
% arbitrary curve framing.

% Usage: [x,y,z]=tubeplot(curve,frame,r,n)

% Constructs a tube around a given curve using a given normal framing.

% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
% Example of use:
% t=linspace(0,2*pi,50);
% curve = [cos(t);sin(t);0.2*(t-pi).^2];
% [~,frame,~] = frenet_frame(p)
% tubeplot(curve,frame,0.1);
% daspect([1,1,1]); camlight;
%
% Arguments:
% curve  [3,N] vector of curve data
% frame  [3,N] vector of normal field data
% r      the radius of the tube
% n      number of points to use on circumference. Defaults to 8

%
% The algorithms fails if you have bends beyond 90 degrees.
% Janus H. Wesenberg, july 2004

  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
    end;
  end;
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end;

  [~,npoints]=size(curve);

  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);
  
  %precalculate cos and sin factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: define the surface according to the given framing
  for k=1:npoints
    convec=cross(frame(:,k),dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
  end;
  
  
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  
  %... and plot:
  if nargout<3, surf(x,y,z); end;
  