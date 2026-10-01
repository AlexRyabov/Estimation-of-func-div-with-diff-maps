function s = NS(n, varargin)
%NS Number to string (short form of num2str).
%
%   s = NS(n)
%   s = NS(n, precision_or_format)

if ~isempty(varargin)
    s = num2str(n, varargin{1});
else
    s = num2str(n);
end
end
