clear;
clc;
close all;

addpath(['..' filesep '..' filesep 'Misc' filesep]); 

linewidth = 1.5;
fontname = 'Times';
fontsize = 24;

figProp = struct( 'size' , fontsize , 'font' ,fontname , 'lineWidth' , linewidth, 'figDim', [1 1 800 600]);
figure

[lambda, RGB] = createSpectrum('1964_full');

image(lambda, 1:size(RGB,1), RGB);

ax = gca;
ax.Visible = 'off';
set(gca,'color','none')

formatFig( gcf ,['.' filesep 'figs' filesep 'spectrum'],'en' , figProp );