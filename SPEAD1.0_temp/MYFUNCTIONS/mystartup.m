%----------------------------------------------------------------------------
% <https://www.mathworks.com/help/matlab/ref/startup.html> 
% <https://www.mathworks.com/help/matlab/ref/matlabrc.html>
% <https://www.mathworks.com/help/matlab/matlab_env/preferences.html> 
%----------------------------------------------------------------------------
% For individual systems USE the startup.m file to customize MATLAB startup. 
% For multiusers systems the matlabrc.m file (located in the matlabroot/toolbox/local folder) is reserved for system admin.
%----------------------------------------------------------------------------
%............................................................................
close all
clear all
format short g
%............................................................................
addpath ~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYFUNCTIONS/
%............................................................................
GenPath = genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYTOOLBOX/');
addpath(GenPath)
%............................................................................
%FIGURE RENDERING FOR MATLAB: 
% <http://matlab.izmiran.ru/help/techdoc/graphics_csh/error_pages/opengl_error_mac.html> 
opengl('neverselect')
set(0,'DefaultFigureRenderer','painters')
set(0,'DefaultFigureToolbar','none')
%............................................................................
%----------------------------------------------------------------------------
% <https://stackoverflow.com/questions/33742132/matlab-figure-rendering-opengl-vs-painters> 
% <http://www.mathworks.com/help/matlab/ref/figure-properties.html#property_renderer> 
% NOTE : Painters is usually faster. OpenGL makes a big difference in the quality of 3-D plots. 
% Especially if you're using lighting or transparency, or if you have a Retina display.
hplot = plot(0,0,'ko','markersize',50,'linewidth',8);
get(gcf,'Renderer') %GCF (get current figure) 
set(gcf,'Renderer')
set(gcf,'Renderer','opengl')
set(gcf,'Renderer','painters')
set(gcf,'Renderer','zbuffer')
set(gcf,'Renderer','none')
%----------------------------------------------------------------------------
%............................................................................
%FIGURE RENDERING FOR OCTAVE: 
% $$$ graphics_toolkit()
% $$$ graphics_toolkit("fltk")
% $$$ graphics_toolkit("gnuplot")
% $$$ available_graphics_toolkits()
%............................................................................
