function f_Set_scales(Scales, varargin)
%F_SET_SCALES Set linear or logarithmic axis scales and the frame style.
%
%   f_Set_scales({'x', 'y'})               linear axes
%   f_Set_scales({'logx', 'logy'})         logarithmic axes
%   f_Set_scales({'logx', 'logy', 'logz'})
%   f_Set_scales(Scales, dir)              tick direction passed to f_FramePlot

if any(strcmp(Scales, 'logx'))
    set(gca, 'XScale', 'log')
end
if any(strcmp(Scales, 'x'))
    set(gca, 'XScale', 'linear')
end
if any(strcmp(Scales, 'logy'))
    set(gca, 'YScale', 'log')
end
if any(strcmp(Scales, 'y'))
    set(gca, 'YScale', 'linear')
end
if any(strcmp(Scales, 'logz'))
    set(gca, 'ZScale', 'log')
end
if any(strcmp(Scales, 'z'))
    set(gca, 'ZScale', 'linear')
end
if isempty(varargin)
    f_FramePlot;
else
    f_FramePlot(varargin);
end
end
