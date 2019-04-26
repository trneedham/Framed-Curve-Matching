function [z1best,w1best,z2best,w2best,gamIbest] = Optimal_Reparam(p1,p2,seeds)

% Input: A pair of closed space curves and some seed positions for p2.
% Output: Optimal choice of p2 which minimizes geod dist, where for each seed position
%         we optimize over rotations and reparameterizations, then take the
%         best choice from all of the seeds. The output is the Stiefel
%         coordinates of the best curve, as well as the optimal
%         reparameterization, gamIbest.

[~,k]=size(seeds);
[~,n]=size(p1);

q1=curve2quat_frenet(p1);
[z1,w1]=quat_to_stiefel(q1);

minDist = 1000;

for j=1:k
    p2j=ShiftF(p2,seeds(j));
    q2j=curve2quat_frenet(p2j);
    [z2j,w2j]=quat_to_stiefel(q2j);
    
    [~,z1new,w1new,z2new,w2new,~,~]=Grassmannian_Geodesic(z1,w1,z2j,w2j);
    
    for j=1:n
        q1j(1,j)=real(z1new(j));
        q1j(2,j)=imag(z1new(j));
        q1j(3,j)=real(w1new(j));
        q1j(4,j)=imag(w1new(j));
    end
    
    for j=1:n
        q2j(1,j)=real(z2new(j));
        q2j(2,j)=imag(z2new(j));
        q2j(3,j)=real(w2new(j));
        q2j(4,j)=imag(w2new(j));
    end
    
    [gam] = DynamicProgrammingQ2(q1j,q2j,0,0);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    p2new = Group_Action_by_Gamma_Coord(p2j,gamI);
    q2new=curve2quat_frenet(p2new);
    [z2new,w2new]=quat_to_stiefel(q2new);
    
    distj=Grassmannian_Geodesic_Dist(z1new,w1new,z2new,w2new);
    
    if distj < minDist
        z1best=z1new;
        w1best=w1new;
        z2best=z2new;
        w2best=w2new;
        gamIbest=gamI;
        minDist=distj;
    end
end


return;


