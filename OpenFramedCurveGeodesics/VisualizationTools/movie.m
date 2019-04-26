% Input is the path geodesic. The output plays a movie of the geodesic.

function m=movie(geod)

k=size(geod,3);

hold off

for i=1:k
    plot3(geod(1,:,i),geod(2,:,i),geod(3,:,i));
    axis equal
    pause(0.1)
end