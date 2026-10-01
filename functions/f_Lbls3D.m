function f_Lbls3D(xLabel, yLabel, varargin)
%F_LBLS3D Set x, y and, optionally, z labels and the title of the current axes.
%
%   f_Lbls3D(xLabel, yLabel)
%   f_Lbls3D(xLabel, yLabel, zLabel)
%   f_Lbls3D(xLabel, yLabel, zLabel, Title)

xlabel(xLabel);
ylabel(yLabel);
if nargin > 2
    zlabel(varargin{1});
end
if nargin > 3
    title(varargin{2});
end
end
