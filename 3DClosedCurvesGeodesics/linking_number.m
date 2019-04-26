function link=linking_number(p1,p2)

[~,n]=size(p1);

dp1=polEdge(p1);
dp2=polEdge(p2);

dt=1/n;

for j=1:n
    for k=1:n
        matrix{j,k}=[dp1(:,j),dp2(:,k),p1(:,j)-p2(:,k)];
    end
end

for j=1:n
    for k=1:n
        difference_norm(j,k)=norm(p1(:,j)-p2(:,k));
    end
end

for j=1:n
    for k=1:n
        integrand(j,k)=det(matrix{j,k})/difference_norm(j,k)^3;
    end
end


link = round((1/(4*pi))*sum(integrand(:))*dt^2);