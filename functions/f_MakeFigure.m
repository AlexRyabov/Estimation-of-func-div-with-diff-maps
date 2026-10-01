function Fg = f_MakeFigure(Fg, Position)
%F_MAKEFIGURE Create (or select) a figure, clear it and set its position.
%
%   Fg = f_MakeFigure(Fg, Position)
%   Fg is a figure number or handle; negative elements of Position keep
%   the current values, e.g. [-1, -1, 800, 600] changes only the size.

Fg = figure(Fg);
clf;
ind = Position < 0;
Position(ind) = Fg.Position(ind);
set(Fg, 'Position', Position);
end
