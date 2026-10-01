function [lData, lFit, mFit, sLegend] = f_PlotAndFit(x, y, FitType, Marker, Color)
%F_PLOTANDFIT Plot data and a fitted trend.
%
%   [lData, lFit, mFit, sLegend] = f_PlotAndFit(x, y, FitType, Marker, Color)
%   plots y vs x (NaN values of y are ignored) and fits
%     FitType = 'lin'   linear model (fitlm)
%     FitType = 'pol3'  third-order polynomial (fitnlm)
%   x can be numeric or datetime (datetime is converted with datenum).
%
%   Outputs: data line, fit line, fitted model and a legend string with the
%   slope and its p-value.

notnan = ~isnan(y);
y = y(notnan);
x = x(notnan);
lData = plot(x, y, Marker, 'Color', Color);
lData.Annotation.LegendInformation.IconDisplayStyle = 'on';
hold on

if isdatetime(x)
    x2fit = datenum(x);
else
    x2fit = x;
end

switch FitType
    case 'lin'
        mFit = fitlm(x2fit, y);
    case 'pol3'
        modelfun = @(b, x) b(1) + b(2) * x + b(3) * x.^2 + b(4) * x.^3;
        mFit = fitnlm(x2fit, y, modelfun, [y(1), 0, 0, 0]);
    otherwise
        error('Unknown fit type: %s', FitType)
end

x2plot = linspace(min(x2fit), max(x2fit), 500)';
lFit = plot(datetime(datestr(x2plot)), mFit.predict(x2plot), 'Color', Color, 'LineWidth', 2);
lFit.Annotation.LegendInformation.IconDisplayStyle = 'off';
sLegend = sprintf('\\alpha=%.2g, p=%.2g', mFit.Coefficients.Estimate(2), mFit.Coefficients.pValue(2));
end
