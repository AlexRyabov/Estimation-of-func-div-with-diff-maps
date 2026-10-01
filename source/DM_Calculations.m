function Result = DM_Calculations(CalcID, ~, Data)
%DM_CALCULATIONS Auxiliary calculations for DM_CompareSpec.
%
%   Result = DM_Calculations('Transp', [], Data)
%   Transposes the species x samples biomass table, converts biomass to
%   relative biomass, adds date, year, month, season, station and its
%   coordinates to each sample, and calculates the Rao functional
%   diversity, species richness and Simpson effective species number.
%
%   Data fields:
%     tPlAbndUns   species x samples biomass table; sample columns are
%                  named <prefix>_<yyyyMMdd>_<station>
%     indAbnd      indexes of the biomass columns in tPlAbndUns
%     tStatCoords  station names and coordinates
%     aEV          diffusion-map eigenvectors (species x i-traits)
%     DiffDist     diffusion distance between species
%
%   Result fields:
%     tPlAbndTranspFreq  samples x species table of relative biomass with
%                        additional columns (see above)
%     indAbndTranspFreq  indexes of the species columns
%
%   The second input (figure number) is not used; it is kept for
%   consistency with DM_Adds_Calculations.

switch CalcID
    case 'Transp'
        tPlAbndUns  = Data.tPlAbndUns;
        indAbnd     = Data.indAbnd;
        tStatCoords = Data.tStatCoords;
        aEV         = Data.aEV;
        DiffDist    = Data.DiffDist;

        % Column names for the transposed table: <ALGAEGROUP>_<SPECIES>
        AlgSpecNames = strcat(tPlAbndUns.ALGAEGROUPS, '_', tPlAbndUns.SPECIES);
        AlgSpecNames = matlab.lang.makeValidName(AlgSpecNames);

        % Extract sampling date and station from the column names
        C = split(tPlAbndUns.Properties.VariableNames(indAbnd), '_');
        C = squeeze(C);
        SampleDate = datetime(C(:, 2), 'InputFormat', 'yyyyMMdd');
        Station = C(:, 3);

        % Station coordinates
        Lat = NaN(size(Station));
        Lon = NaN(size(Station));
        for ist = 1:size(tStatCoords, 1)
            tmp_ind = strcmp(char(tStatCoords.StationNames(ist)), Station);
            if any(tmp_ind)
                Lat(tmp_ind) = tStatCoords.WGS_84_Y_DecimalDegree(ist);
                Lon(tmp_ind) = tStatCoords.WGS_84_X_DecimalDegree(ist);
            end
        end

        % Transpose: samples in rows, species in columns
        Xc = table2cell(tPlAbndUns(:, indAbnd));
        tPlAbndTranspFreq = cell2table(Xc', 'RowNames', tPlAbndUns.Properties.VariableNames(indAbnd), ...
            'VariableNames', AlgSpecNames);
        indAbndTranspFreq = 1:size(tPlAbndTranspFreq, 2);

        % Total biomass and relative biomass of each species
        tPlAbndTranspFreq.Total = sum(tPlAbndTranspFreq{:, indAbndTranspFreq}, 2);
        tPlAbndTranspFreq{:, indAbndTranspFreq} = tPlAbndTranspFreq{:, indAbndTranspFreq} ./ ...
            repmat(tPlAbndTranspFreq.Total(:), 1, length(indAbndTranspFreq));

        % Date, year, month, coordinates, season and station
        tPlAbndTranspFreq.Date  = SampleDate;
        tPlAbndTranspFreq.Year  = year(SampleDate);
        tPlAbndTranspFreq.Month = month(SampleDate);
        tPlAbndTranspFreq.Lat   = Lat;
        tPlAbndTranspFreq.Lon   = Lon;
        Seasons = {'Winter', 'Winter', 'Spring', 'Spring', 'Spring', 'Summer', ...
                   'Summer', 'Summer', 'Autumn', 'Autumn', 'Autumn', 'Winter'};
        tPlAbndTranspFreq.Season  = Seasons(tPlAbndTranspFreq.Month)';
        tPlAbndTranspFreq.Station = Station;

        % Diversity indexes of each sample
        X = tPlAbndTranspFreq{:, indAbndTranspFreq};
        tPlAbndTranspFreq.FuncDivRao = f_DiversityMetrics(X, aEV, 'RaoDist', DiffDist, 2)';
        tPlAbndTranspFreq.Richness   = f_DiversityMetrics(X, [], 'Richness', [], 2)';
        tPlAbndTranspFreq.SimpsonESN = f_DiversityMetrics(X, [], 'SimpsonESN', [], 2)';

        Result.tPlAbndTranspFreq = tPlAbndTranspFreq;
        Result.indAbndTranspFreq = indAbndTranspFreq;
    otherwise
        error('CalcID %s is not defined in DM_Calculations', CalcID)
end
end
