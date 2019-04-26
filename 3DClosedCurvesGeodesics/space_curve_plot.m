% This gives a plot of a space curve as usual, just with less typing.
% Input is a space curve, entered as a 3xn matrix, output is the plot.
% Options can be adjusted here per user preferences.

function [curve]=space_curve_plot(p)

curve=plot3(p(1,:),p(2,:),p(3,:),'LineWidth',3);
axis equal