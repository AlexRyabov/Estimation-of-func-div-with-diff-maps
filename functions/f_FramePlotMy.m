function f_FramePlotMy(varargin)
%F_FRAMEPLOTMY Frame style used in the figures of the paper.
%
%   f_FramePlotMy             box, minor ticks, ticks inside
%   f_FramePlotMy(dir)        tick direction 'in' or 'out'
%   f_FramePlotMy(dir, len)   tick direction and tick length

f_FramePlot('in')
box off
switch nargin
    case 0
        set(gca, 'YMinorTick', 'on', 'box', 'on', 'XMinorTick', 'on', 'TickDir', 'in', 'Layer', 'top', 'TickLength', [0.02, 0.02]);
    case 1
        set(gca, 'YMinorTick', 'on', 'XMinorTick', 'on', 'TickDir', char(varargin{1}), 'Layer', 'top', 'TickLength', [0.02, 0.02]);
    case 2
        set(gca, 'YMinorTick', 'on', 'XMinorTick', 'on', 'TickDir', char(varargin{1}), 'Layer', 'top', 'TickLength', [varargin{2}, varargin{2}]);
end
end
