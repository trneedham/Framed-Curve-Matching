% Input is the path geodesic. The output plays a movie of the geodesic.

function m=movie3D(geod)

k=size(geod,3);

hold off

for i=1:k
    space_curve_plot(geod(:,:,i));
    pause(0.1)
end