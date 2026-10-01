function [Date, Station] = f_GetStationDate(s)
%F_GETSTATIONDATE Sampling date and station from a sample name.
%
%   [Date, Station] = f_GetStationDate('C_19930514_6') returns
%   Date = 14-May-1993 and Station = '6'. The sample names are the column
%   names of the biomass table in DM_CompareSpec. For a cell array of names
%   a datetime vector and a cell array of stations are returned.

if iscell(s)
    Date = NaT(length(s), 1);
    Station = cell(length(s), 1);
    for i = 1:length(s)
        [Date(i), Station{i}] = f_GetStationDate(s{i});
    end
else
    S = strsplit(s, '_');
    Date = datetime(S{2}, 'InputFormat', 'yyyyMMdd');
    Station = S{3};
end
end
