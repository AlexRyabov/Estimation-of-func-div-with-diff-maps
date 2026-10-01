function f_FramePlot(varargin)
%F_FRAMEPLOT Box, minor ticks and tick style of the current axes.
%
%   f_FramePlot             ticks inside
%   f_FramePlot(dir)        tick direction 'in' or 'out'
%   f_FramePlot(dir, len)   tick direction and tick length

box on
switch nargin
    case 0
        set(gca, 'YMinorTick', 'on', 'XMinorTick', 'on', 'TickDir', 'in', 'Layer', 'top', 'TickLength', [0.02, 0.02]);
    case 1
        set(gca, 'YMinorTick', 'on', 'XMinorTick', 'on', 'TickDir', char(varargin{1}), 'Layer', 'top', 'TickLength', [0.02, 0.02]);
    case 2
        set(gca, 'YMinorTick', 'on', 'XMinorTick', 'on', 'TickDir', char(varargin{1}), 'Layer', 'top', 'TickLength', [varargin{2}, varargin{2}]);
end
end
