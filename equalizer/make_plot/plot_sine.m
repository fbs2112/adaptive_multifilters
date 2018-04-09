clear;
clc;
close all

t = 0:0.01:5;

plot(sin(2*pi*t),'k','linewidth',2)

ax = gca;
ax.Visible = 'off';

print('sine','-dpng');