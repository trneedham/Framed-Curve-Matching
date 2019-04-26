function qn = Group_Action_by_Gamma_Coord_Quaternions(q,gamma)

[n,T] = size(q);

for j=1:n
    qn(j,:) = spline(linspace(0,1,T) , q(j,:) ,gamma);
end

dgamma = gradient(gamma);

for j=1:T
    qn(:,j) = sqrt(dgamma(j))*qn(:,j);
end

return;