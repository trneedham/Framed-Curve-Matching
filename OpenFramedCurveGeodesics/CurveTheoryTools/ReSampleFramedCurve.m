% This takes a curve and frame and resamples them uniformly N times. 
% Modified from the code of Srivastava et al.

function [Xn,Vn] = ReSampleFramedCurve(X,V,N)

    [n,T] = size(X);
    del(1) = 0;
    for r = 2:T
        del(r) = norm(X(:,r) - X(:,r-1));
    end
    cumdel = cumsum(del)/sum(del);   
    
    newdel = [0:N-1]/(N-1);
    
    for j=1:n
        Xn(j,:) = interp1(cumdel,X(j,1:T),newdel,'linear');
        Vn(j,:) = interp1(cumdel,V(j,1:T),newdel,'linear');
    end
    
    for j=1:N
        if norm(Vn(:,j)) == 0
            Vn(:,j)=Vn(:,j);
        else
            Vn(:,j)=Vn(:,j)/norm(Vn(:,j));
        end
    end