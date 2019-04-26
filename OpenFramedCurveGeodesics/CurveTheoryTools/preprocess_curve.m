% Takes a curve p and preprocesses it by: 
% -shifting to begin at the origin
% -scaling so that it has total length 1.

function pnew=preprocess_curve(p)

[d,n]=size(p);

% Shift to begin at the origin.
for i=1:d
   pnew2(i,:)=p(i,:)-p(i,1);
end

% Now we calculate the length of the curve. 

curveLength=arclength(pnew2(1,:),pnew2(2,:),pnew2(3,:));

% Rescale so that the curve has total length 1.

pnew=pnew2/curveLength;
    
    