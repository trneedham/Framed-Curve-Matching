% This function takes a space curve given by a list of points parameterized 
% by time and outputs the discrete Frenet frame --- in particular, the normal
% and binormal vectors. V is normal vector and W is the binormal.

% This algorithm is modified to handle vanishing curvature by a 'brute
% force' approach. Vanishing curvature manifests as consective edges pointing in 
% the same direction---when this occurs we simply parallel transport the frames.
% This requires the curve to have nonvanishing curvature somewhere to start with.
% Then the rest of the frame is defined based on the first nonvanishing curvature
% point. 

function [U,V,W] = frenet_frame(p)

[~,n] = size(p);

% Discrete tangent curve
edges=polEdge(p);

% Discrete unit tangents
for i=1:n
    U(:,i)=edges(:,i)/norm(edges(:,i));
end

% Define a list containing all cross products.
for i=1:n-1
  CrossProds(:,i)=cross(edges(:,i),edges(:,i+1));
end

% Now we search for a point of nonvanishing curvature.
for i=1:n-1
  CrossProdsNorms(i)=norm(CrossProds(:,i));
end

f=find(CrossProdsNorms>0.000001,1);

% If f is empty, then p is a straight line and we define the frame by
% choosing a random vector and projecting to the normal plane.
if isempty(f)==1
    rand_vect=rand(3,1);
    tang_vect=p(:,2)-p(:,1);
    proj_vect=dot(rand_vect,tang_vect)/norm(tang_vect)^2*tang_vect;
    norm_vect=(rand_vect-proj_vect)/norm(rand_vect-proj_vect);
    
    for i=1:n
        normVect(:,i)=norm_vect;
        binorm(:,i)=cross(tang_vect,norm_vect);
    end
else
    % Now we define the binormals via the cross products, or parallel transporting 
    % if the cross product is zero. We start at the nonvanishing curvature point 
    % and work outwards.
    binorm(:,f)=CrossProds(:,f)/CrossProdsNorms(f);

    for i=f+1:n-1
        if (CrossProdsNorms(i)>0.000001)
            binorm(:,i)=CrossProds(:,i)/CrossProdsNorms(i);
        else
            binorm(:,i)=binorm(:,i-1);
        end
    end

    % Extend the binormal list so that it has the same number of entries as edges
    binorm(:,n)=binorm(:,n-1);

    for i=1:f-1
        if (CrossProdsNorms(f-i)>0.000001)
            binorm(:,f-i)=CrossProds(:,i)/CrossProdsNorms(i);
        else
            binorm(:,f-i)=binorm(:,f-i+1);
        end
    end
    
    % Finally, we complete the frame by defining the normal vectors via cross product.
    for i=1:n
        normVect(:,i)=cross(binorm(:,i),edges(:,i)/norm(edges(:,i)));
    end
end
  
V=normVect;
W=binorm;
