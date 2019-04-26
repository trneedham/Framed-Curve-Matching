function dist = RMS_dist(p1,p2)

p1=preprocess_curve(ReSampleCurve(p1,50));
p2=preprocess_curve(ReSampleCurve(p1,50));

dist=1000;

for j=1:50
    p2n = ShiftF(p2,j-1);
    A = p1*p2n';
        [U,~,V] = svd(A);
        if det(A)> 0
            Ot = U*V';
        else
            if (size(p1,1)==2)
                Ot = U*([V(:,1) -V(:,2)])';
            else
                Ot = U*([V(:,1) V(:,2) -V(:,3)])';
            end
        end
        p2n = Ot*p2n;
        
        distn=norm(p1-p2n);
        
        if distn<dist;
            dist=distn;
        end
end