close all
clc
clear 

figure
h=patch([0 1 0 1], [0 1 1 0], 'r');
x=.5*ones(10,1);
y=linspace(0,1, 10);
hold on
hp=plot(x,y);

%Returns handles to the patch and line objects
chi=get(gca, 'Children');
%Reverse the stacking order so that the patch overlays the line
set(gca, 'Children',flipud(chi))