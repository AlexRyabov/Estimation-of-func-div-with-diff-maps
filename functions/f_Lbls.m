function f_Lbls(xLabel, yLabel, varargin)
%F_LBLS Set x and y labels and, optionally, the title of the current axes.
%
%   f_Lbls(xLabel, yLabel)
%   f_Lbls(xLabel, yLabel, Title)

xlabel(xLabel);
ylabel(yLabel);
if ~isempty(varargin)
    title(varargin{1})
end
end
