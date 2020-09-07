%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: myLatexYticks.m
%
%Usage: [hy] = myLatexYticks(h,ticky,tickposy,roty,offset,varargin);
%
%Description: Replace or appends XTickLabels and YTickLabels of axis handle
%             h with input tickx and ticky array
%
%***NOTE!***: BE SURE TO DELETE ANY PREVIOUS TEXT OBJECTS CREATED BY THIS
%             FUNCTION BEFORE RUNNING THIS ON THE SAME FIGURE TWICE
%
%Required Inputs:
%       h        : handle of axis to change tick labels (can use gca)
%
%Optional Inputs
%       ticky    : cell array of tick labels or string to append to current
%                  labels (Can use [] or not specify to ignore) 
%       tickposy : Vector of y positions where you want the tick labels
%                  (Can use [] or not specify to ignore) 
%       roty     : Number of degrees to rotate y tick labels
%                  (Can use [] or not specify to ignore) Default = 0.0
%       offset   : Label offsets from axis in fraction of total range
%                  (Can use [] or not specify to ignore) Default = 0.0
%
%Optional Inputs:%                
%                Any standard text formatting parameters such as 
%                'FontSize','FontWeight',etc.
%                Use the same way you would in a set command after putting
%                in the required input values.
%
%Outputs:
%        hy: handle of text objects created for YTickLabels
%
%Function Calls:
%               None
%
%Required Data Files:
%                    None
%
%Change Log:
%           08/19/2007: Origin Version Created by Alex Hayes
%                       (hayes@gps.caltech.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN FUNCTION ("Function FORMAT_TICKS.m"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hy] = myLatexXticks(h,ticky,tickposy,roty,offset,varargin)

%define axis text offset (percentage of total range)
if ~exist('offset','var');
    offset = 0.02;
elseif length(offset) == 0;
    offset = 0.02;
end;

%make sure the axis handle input really exists
if ~exist('h','var');
    h = gca;
    warning(['Axis handle NOT Input, Defaulting to Current Axes num2str(h)']);
elseif length(h) == 0;
    h = gca;
    warning(['Axis Handle NOT Input, Defaulting to Current Axes num2str(h)']);
elseif ~ishandle(h(1))
    warning(['Input (' num2str(h(1)) ') is NOT an axis handle, defaulting to current axis, num2str(h)']);
    h = gca;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BEGIN: NOW THE Y-AXIS TICK LABELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%only move forward if we are doing anything to the yticks
if ~exist('ticky');
    hy = -1;
elseif length(ticky)==0;
    hy = -1;
else;
    %fix the YTickLabels if they have been erased in the past
    if length(get(h,'YTickLabel'))==0;
        set(h,'YTickLabel',get(h,'YTick'));
    end;
    %set the ytick positions if entered
    if exist('tickposy','var');
        if length(tickposy) > 0;
            set(h,'YTick',tickposy);
            set(h,'YTickLabel',tickposy);
        end;
    end;
    %make sure the xtick positions are in the xlimit range
    if exist('tickposy','var');
        if length(tickposy) > 0;
            lim = get(h,'YLim');
            if lim(1) > min(tickposy);
                lim(1) = min(tickposy);
            end;
            if lim(2) < max(tickposy);
                lim(2) = max(tickposy);
            end;
            set(h,'YLim',lim);
        end;
    end;
    %get the tick labels and positions if the user did not input them
    if ~exist('ticky','var');
        ticky = get(h,'YTickLabel');
        if ischar(ticky);
            temp = ticky;
            ticky = cell(1,size(temp,1));
            for j=1:size(temp,1);
                ticky{j} = strtrim( temp(j,:) );
            end;
        end;
        append = '^{\circ}';
        for j=1:length(ticky);
            ticky{j} = [ticky{j} append];
        end;
    elseif length(ticky) == 0;
        ticky = get(h,'YTickLabel');
        if ischar(ticky);
            temp = ticky;
            ticky = cell(1,size(temp,1));
            for j=1:size(temp,1);
                ticky{j} = strtrim( temp(j,:) );
            end;
        end;
        append = '^{\circ}';
        for j=1:length(ticky);
            ticky{j} = [ticky{j} append];
        end;
    elseif isstr(ticky);
        append = ticky;
        ticky = get(h,'YTickLabel');
        if ischar(ticky);
            temp = ticky;
            ticky = cell(1,size(temp,1));
            for j=1:size(temp,1);
                ticky{j} = strtrim( temp(j,:) );
            end;
        end;
        if strcmp(append(1),'$');
            for j=1:length(ticky);
                ticky{j} = ['$' ticky{j} append(2:end)];
            end;
        else;
            for j=1:length(ticky);
                ticky{j} = [ticky{j} append];
            end;
        end;
    elseif ~iscell(ticky );
        warning(['Input TICKY variable is not a compatible string ' ...
            'or cell array! Returning...']);
        return;
    end;
    %find out if we have to use the LaTex interpreter
    temp = ticky{1};
    if strcmp(temp(1),'$');
        latex_on = 1;
    else;
        latex_on = 0;
    end;
    %erase the current tick label
    set(h,'YTickLabel',{});
    %get the x tick positions if the user did not input them
    if ~exist('tickposy','var');
        tickposy = get(h,'YTick');
    elseif length(ticky) == 0;
        tickposy = get(h,'YTick');
    end;
    %get the x tick positions if the user did not input them
    if ~exist('tickposx','var');
        tickposx = get(h,'YTick');
    elseif length(tickposx) == 0;
        tickposx = get(h,'XTick');
    end;
    %set the new tick positions
    set(h,'YTick',tickposy);
%    set(h,'XTick',tickposx);
    %check the lengths of the xtick positions and xtick labels
    l1 = length(ticky);
    l2 = length(tickposy);
    if l1==0;
        set(h,'YTickLabel',ticky);
    end;
    if l1~=l2;
        disp(['Length of YTick = ' num2str(length(tickposy))]);
        disp(['Length of YTickLabel = ' num2str(length(ticky))]);
        if l2 < l1;
            warning(['Reducing Length of YTickLabel!']);
        else;
            warning(['Reducing Length of YTick!']);
        end;
        l3 = min([l1,l2]);
        ticky = ticky{1:l3};
        tickposy = tickposy(1:l3);
    end;
    %set rotation to 0 if not input
    if ~exist('roty','var');
        roty = 0;
    elseif length(roty) == 0;
        roty = 0;
    end;
    %Convert the cell labels to a character string
    %ticky = char(ticky);
    ticky = cellstr(ticky);
    %Make the YTICKS!
    lim = get(h,'XLim');
    if min(tickposx) < lim(1);
        lim(1) = min(tickposx);
    end;
    if max(tickposx) > lim(2);
        lim(2) = max(tickposx);
    end;
    if roty == 0;
        if latex_on;
            hy = text(...
                repmat(lim(1)-offset*(lim(2)-lim(1)),length(tickposy),1),...
                tickposy,...
                ticky,'VerticalAlignment','middle',...
                'HorizontalAlignment','right','rotation',roty,'interpreter','LaTex');
        else;
            hy = text(...
                repmat(lim(1)-offset*(lim(2)-lim(1)),length(tickposy),1),...
                tickposy,...
                ticky,'VerticalAlignment','middle',...
                'HorizontalAlignment','right','rotation',roty);
        end;
    elseif roty < 180;
         if latex_on;
            hy = text(...
                repmat(lim(1)-offset*(lim(2)-lim(1)),length(tickposy),1),...
                tickposy,...
                ticky,'VerticalAlignment','middle',...
                'HorizontalAlignment','right','rotation',roty,'interpreter','LaTex');
        else;
            hy = text(...
                repmat(lim(1)-offset*(lim(2)-lim(1)),length(tickposy),1),...
                tickposy,...
                ticky,'VerticalAlignment','middle',...
                'HorizontalAlignment','right','rotation',roty);
        end;
    else;
          if latex_on;
            hy = text(...
                repmat(lim(1)-offset*(lim(2)-lim(1)),length(tickposy),1),...
                tickposy,...
                ticky,'VerticalAlignment','middle',...
                'HorizontalAlignment','right','rotation',roty,'interpreter','LaTex');
        else;
            hy = text(...
                repmat(lim(1)-offset*(lim(2)-lim(1)),length(tickposy),1),...
                tickposy,...
                ticky,'VerticalAlignment','middle',...
                'HorizontalAlignment','right','rotation',roty);
        end;
    end;
    %Get and set the text size and weight
    set(hy,'FontSize',get(h,'FontSize'));
    set(hy,'FontWeight',get(h,'FontWeight'));

    %Set the additional parameters if they were input
    if length(varargin) > 2;
        command_string = ['set(hy'];
        for j=1:2:length(varargin);
            command_string = [command_string ',' ...
                '''' varargin{j} ''',varargin{' num2str(j+1) '}'];
        end;
        command_string = [command_string ');'];
        eval(command_string);
    end;
end;
        
