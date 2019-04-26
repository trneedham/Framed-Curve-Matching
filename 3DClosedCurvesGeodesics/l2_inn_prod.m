function z=l2_inn_prod(z1,z2)

% Input is a pair of complex vectors, treated as complex-valued paths.
% Output is the complex L^2 inner product of the complex paths.

[~,n]=size(z1);
dt=1/n;


for j=1:n
    inn_prods(j)=z1(j)*conj(z2(j));
end

z=sum(inn_prods)*dt;