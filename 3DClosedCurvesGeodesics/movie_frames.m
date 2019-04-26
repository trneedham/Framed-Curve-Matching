% Input is a geodesic of a framed path. geod is the geodesic of the paths 
% and pushoff_geod is a geodesic of the framing. Output plays a movie of
% the geodesic.

function m=movie_frames(geod,pushoff_geod)

k=size(geod,3);

clf

for i=1:k-1
    hold on
    view(45,45)
    space_curve_plot(geod(:,:,i));
    space_curve_plot(pushoff_geod(:,:,i));
    pause(0.1)
    hold off
    clf
end

hold on 
view(45,45)
space_curve_plot(geod(:,:,i));
space_curve_plot(pushoff_geod(:,:,i));
hold off