function [zn,wn] = Group_Action_by_Gamma_Coord_Stiefel(z,w,gamma)

[~,T] = size(z);

zn(:) = spline(linspace(0,1,T) , z(:) ,gamma);
wn(:) = spline(linspace(0,1,T) , w(:) ,gamma);  

dgamma = gradient(gamma);

for j=1:T
    zn(j) = sqrt(dgamma(j))*zn(j);
    wn(j) = sqrt(dgamma(j))*wn(j);
end

zn=zn/norm(zn);
wn=wn/norm(wn);

return;