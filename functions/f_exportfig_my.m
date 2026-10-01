function f_exportfig_my(fg, filename, part)
%F_EXPORTFIG_MY Export a figure with a resolution scaled to its width.
%
%   f_exportfig_my(fg, filename, part)
%   part is the fraction of the page width occupied by the figure:
%   part = 1 gives an image about 2500 pixels wide, 0.5 about 1250 pixels.

w = fg.Position(3);
r = ceil(2500 / (w * 0.0093) * part);
exportgraphics(fg, filename, 'Resolution', r);
end
