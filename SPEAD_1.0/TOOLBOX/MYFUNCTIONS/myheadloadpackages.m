function [mypackages] = myheadloadpackages(doPlotting)

%****************************************************************************
% PROGRAM: MYHEADLOADPACKAGES.M
%
% ## Use: myheadloadpackages
% Use: [subplot_funhan,colorbar_funhan,verticales,horizontal] = myheadloadpackages(doPlotting);
%
%****************************************************************************

%%%%%%%%%%%%%%%%%%
%GLOBAL VARIABLES:
%%%%%%%%%%%%%%%%%%
%============================================================================
% (THIS MUST BE LOCATED * BEFORE * DECLARING THE FUNCTIONS) 
%----------------------------------------------------------------------------
% $$$ global subplot_funhan
% $$$ global colorbar_funhan verticales horizontal 
%----------------------------------------------------------------------------
%============================================================================

%%%%%%%%%%%%%%%
%LOAD PACKAGES:
%%%%%%%%%%%%%%%
%============================================================================
% LOAD COMMON PACKAGES:
%----------------------------------------------------------------------------
%more off
%close all
%clear all
%format short g
%----------------------------------------------------------------------------
% $$$ addpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYFUNCTIONS/');
%----------------------------------------------------------------------------
% $$$ addpath(genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYTOOLBOX/'));
% $$$ addpath(genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MYTOOLBOX/MATLAB-LIBRERY/'));
%----------------------------------------------------------------------------
% $$$ %%addpath(genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MEXTOOLBOX/NETCDF-MATLAB6p5/'));
% $$$ %%addpath(genpath('~/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/MEXTOOLBOX/NETCDF-MATLAB2011a-64bit/'));
%----------------------------------------------------------------------------
%============================================================================
%----------------------------------------------------------------------------
% LOAD MATLAB VS. OCTAVE PACKAGES:
% <https://wiki.octave.org/Statistics_package> %Missing Functions of Matlab in Octave.
%----------------------------------------------------------------------------
% # 1ST WAY:
[retval] = myisoctave;
if     retval == 0 %MATLAB
    disp('** This is Matlab **') 
elseif retval == 1 %OCTAVE
    disp('** This is Octave **') 
end
%----------------------------------------------------------------------------
% # 2ND WAY (NICER):
[uiIsOctave,uiIsMatlab] = myisoctavematlab
if     uiIsMatlab == 1 %MATLAB
    disp('** This is MATLAB **') 
    %%FunctionName_Colorbar = 'colorbar';
    FunctionName_Colorbar = 'mycolorbarMatlab6p5';
    colorbar_location = get(eval(FunctionName_Colorbar),'Orientation');
    colorbar_location_vertical = 'vertic';
    colorbar_location_horizont = 'horizontal';
    close all
    FunctionHandle_Colorbar = str2func(FunctionName_Colorbar);
    colorbar_funhan = FunctionHandle_Colorbar;
    opengl('neverselect')
    %%opengl('hardware')
    %%opengl('software')
    feature('accel','off') %JIT accelarator -- disabled (off) or enabled (on)
    jit_enable_asking = feature('accel') %%JIT accelarator when disabled (0) or enabled (1)
elseif uiIsOctave == 1 %OCTAVE
    disp('** This is OCTAVE **') 
    %----------------------------------------------------------------------------
    % OCTAVE NOTE_1:
    % >> pkg install -forge odepkg %INSTALL DE ODE PACKAGE
    % >> pkg("load","odepkg") %LOAD THE ODE PACKAGE (1st way)
    % >> pkg load odepkg %LOAD THE ODE PACKAGE (2nd way) 
    %----------------------------------------------------------------------------
    % OCTAVE NOTE_2:
    % $ emacs /etc/octave.conf
    % > pkg("load","odepkg") %LOAD THE ODE PACKAGE EVERY TIME OCTAVE IS CALLED.
    %----------------------------------------------------------------------------
    %% pkg('load','odepkg'); %load package "ode" (for ode solving) -- NOT NEED IT AFTER OCTAVE 4.4.1 (2018)
    pkg('load','optim') %load package "optim" (for curve fitting)
    pkg('load','financial') %load package "financial" (for nanmean)
    pkg('load','statistics'); %load package "statistics" (for nmds)
    % <https://stackoverflow.com/questions/41040999/speed-up-printing-of-mesh-in-octave>
    %%FunctionName_Colorbar = 'colorbar';
    FunctionName_Colorbar = 'mycolorbarOctave3p4';
    if doPlotting
        colorbar_location = get(eval(FunctionName_Colorbar),'Location');
    end
    colorbar_location_vertical = 'EastOutside';
    colorbar_location_horizont = 'SouthOutside';
    close all
    FunctionHandle_Colorbar = str2func(FunctionName_Colorbar);
    colorbar_funhan = FunctionHandle_Colorbar;
    available_graphics_toolkits()
    if doPlotting
        graphics_toolkit('fltk') %FASTER FOR IMAGESC BUT SLOWER FOR PLOTTING
        %%graphics_toolkit('gnuplot') %FASTER FOR PLOT BUT SLOWER FOR IMAGESC
    end
    gnuplot_binary('gnuplot')
    %jit_enable(0)  %JIT accelarator -- disabled (0) or enabled (1)
    %%jit_enable_asking = jit_enable %%JIT accelarator when disabled (0) or enabled (1)
end
%============================================================================
%MYCOLORBAR ORIENTATION:
%----------------------------------------------------------------------------
% $$$ verticales = 'EastOutside';
% $$$ horizontal = 'SouthOutside';
verticales = colorbar_location_vertical;
horizontal = colorbar_location_horizont;
%----------------------------------------------------------------------------
%MYSUBPLOT FUNCTION HANDLE:
%%FunctionName_subplot = 'subplot';
FunctionName_subplot = 'mysubplotmit';
FunctionHandle_subplot = str2func(FunctionName_subplot);
subplot_funhan = FunctionHandle_subplot;
%----------------------------------------------------------------------------
%JUST CHECKING IF PLOTTING (SUBPLOT AND COLORBAR) IS WORKING OKAY:
if doPlotting
    figure(1)
    subplot_funhan(2,2,1);
    hcbar = colorbar_funhan(verticales); 
    subplot_funhan(2,2,2);
    hcbar = colorbar_funhan(horizontal);
    %----------------------------------------------------------------------------
    subplot_funhan
    colorbar_funhan
    disp('** Just checking MATLAB / OCTAVE plotting -- please wait 1 sec **')
    pause(0.5)
    close all
end
%----------------------------------------------------------------------------
%============================================================================
%OUTPUTS:
%----------------------------------------------------------------------------
mypackages.subplot  = subplot_funhan;
mypackages.colorbar = colorbar_funhan;
mypackages.verticales = verticales;
mypackages.horizontal = horizontal;
%----------------------------------------------------------------------------
%============================================================================
return

