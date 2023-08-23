% f_figureplace(figHandles,m,n) OR
% f_figureplace(figHandles,m,n,screenNumber)
%
% This function places figures in m by n grid on the screen. If there are
% multiple screens, and the optional argument screenNumber is not
% specified, the figures are placed on the primary screen. If it is
% specified, then the figures are placed on the specified screen.
%
% fig: A vector of figure handles
% m: number of rows of the grid
% n: number of columns of the grid

function f_figureplace(figHandles,m,n,varargin)
    if nargin > 4
        error("The function accepts at most 4 input arguments");
    elseif nargin > 3
        scnNum=varargin{1};
    else
        scnNum=1;
    end
    r=groot;
    numScreens=size(r.MonitorPositions,1);
    if numScreens<scnNum
        error("There are only" + numScreens +...
            "screens, cannot place figures on screen" + scnNum);
    end
    scnWidthOffset = r.MonitorPositions(scnNum,1)-1;
    scnHeightOffset = r.MonitorPositions(scnNum,2)-1;
    scnWidth = r.MonitorPositions(scnNum,3);
    scnHeight = r.MonitorPositions(scnNum,4);
    if length(figHandles)>m*n
        error('%d figures can not be fitted in %d x %d grid',length(figHandles),m,n);
    end
    for i = 1:m
        for j = 1:n
            if (j+n*(i-1)) <= length(figHandles)
                figPos=[round((j-1)*scnWidth/n) round((m-i)*scnHeight/m) round(scnWidth/n) round(scnHeight/m)] + [scnWidthOffset scnHeightOffset 0 0];
                set(figHandles((j+n*(i-1))),'OuterPosition',figPos);
            end
        end
    end
end