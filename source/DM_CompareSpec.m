%% DM_CompareSpec - functional diversity of Baltic Sea phytoplankton
% Diffusion-map analysis of the Lithuanian phytoplankton monitoring data
% (Baltic Sea, 1993-2015) as presented in
%
%   Ryabov A. et al. (2022) Estimation of functional diversity and species
%   traits from ecological monitoring data. PNAS 119(43), e2118156119.
%   https://doi.org/10.1073/pnas.2118156119
%
% Workflow
%   1. Load phytoplankton biomass and environmental data (../data).
%   2. Select stations and aggregate biomass per species and sample.
%   3. Assign environmental conditions (+-7 days) to each sample.
%   4. Compute the biomass-weighted mean environment of each species.
%   5. Diffusion map: species similarity -> eigenvectors (i-traits).
%   6. Functional diversity (Rao), richness and Simpson ESN per sample.
%   7. Figures (main text and SI).
%
% Main variables
%   tPlAbnd        long table: one row per species, station and date
%   tPlAbndUns     species (rows) x samples (columns) biomass table;
%                  sample columns are named <prefix>_<yyyyMMdd>_<station>
%   tPlAbndUnsEnv  tPlAbndUns + biomass-weighted mean environment
%   ev             eigenvalues of the diffusion map
%   aEV            eigenvectors scaled by 1/eigenvalue (species x i-traits)
%   aS             sparsified species similarity matrix
%   DiffDist       diffusion distance between species
%   tPlAbndTranspFreq  samples (rows) x species relative biomass, plus
%                  date, station, functional diversity, richness, ESN
%
% Run this script from any folder; paths are resolved relative to it.

clear; close all

%% Parameters
% Station selection:
%   '60vis'  - stations visited at least 60 times (used in the paper)
%   '20vis'  - stations visited at least 20 times
%   'Regs'   - regularly visited stations
%   '2004on' - stations sampled after 2004 (broader monitoring)
%   'All'    - all stations
DataID = '60vis';

% Optional additional plots, e.g. {'VisStations'} shows the station visits
PlotsID = {};

% Save selected figures as PNG files into ../figs
SaveFigures = true;

%% Paths
srcDir = fileparts(mfilename('fullpath'));
if isempty(srcDir)          % e.g. when the script is run section by section
    srcDir = pwd;
end
rootDir = fileparts(srcDir);
dataDir = fullfile(rootDir, 'data');
figDir  = fullfile(rootDir, 'figs');
addpath(srcDir, fullfile(rootDir, 'functions'));
if SaveFigures && ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Load data
load(fullfile(dataDir, 'LT_phytoplankton_Cleaned_Plankton.mat'), '-mat');   % tPlAbnd
load(fullfile(dataDir, 'Chem_Hidrology_TempSal_Avrg_1-10.mat'), '-mat');    % tTempSalDataAvgd
load(fullfile(dataDir, 'Chem_Hidrology_Chem_Avrg_1-10m.mat'), '-mat');      % tChemDataAvgd

tStatCoords = readtable(fullfile(dataDir, 'LT_phytoplankton_Cleaned_StationCoordinates.xlsx'));
tStatCoords.StationNames = categorical(tStatCoords.StationNames);

iFgNr = 50;  % number of the next figure created by DM_Adds_Calculations

%% Select stations
switch DataID
    case 'Regs'
        lStations = {'1B', '2', '20', '3', '4', '5', '6', '64', '7'};
    case '60vis'
        lStations = {'1B', '2', '20', '3', '4', '5', '6', '64', '65', '7'};
    case '20vis'
        lStations = {'1', '1B', '2', '20', '20A', '3', '37', '4', '46', '5', '6', '64', ...
            '65', '6B', '7', 'B-1', 'B-2', 'B-3', 'B-4', 'B-5', 'B-6', 'B-7', 'B-8', 'N-6'};
    case '2004on'
        lStations = {'1B', '2', '20', '3', '4', '5', '6', '64', '7', 'B-1', 'B-4', 'N-6'};
    case 'All'
        lStations = {};
    otherwise
        error('Unknown DataID: %s', DataID)
end

if ~isempty(lStations)
    if strcmp(DataID, '2004on')
        dtStart = datetime(2004, 1, 1);
    else
        dtStart = [];
    end
    tPlAbnd          = SelectStations(tPlAbnd, lStations, dtStart);
    tChemDataAvgd    = SelectStations(tChemDataAvgd, lStations, dtStart);
    tTempSalDataAvgd = SelectStations(tTempSalDataAvgd, lStations, dtStart);
end

tPlAbnd.STATION = categorical(upper(tPlAbnd.STATION));
tPlAbnd.AREA    = categorical(upper(tPlAbnd.AREA));
tPlAbnd.SEASON  = categorical(upper(tPlAbnd.SEASON));
tPlAbnd.Abnd    = tPlAbnd.SPECIESBIOMASS__G_L;   % biomass is used as abundance

tChemDataAvgd.STATION    = categorical(tChemDataAvgd.STATION);
tTempSalDataAvgd.STATION = categorical(tTempSalDataAvgd.STATION);

% Remove unrealistic salinity values
tTempSalDataAvgd.SALINITY__(tTempSalDataAvgd.SALINITY__ > 7.9) = NaN;

%% Aggregate biomass per species and sample
% Set zero where biomass is not defined (4 records)
tPlAbnd.Abnd(isnan(tPlAbnd.Abnd)) = 0;

% Total biomass of each species in each sample
GroupVars = {'STATION', 'Date', 'DateStat', 'AREA', 'SEASON', 'DEPTH_M', 'ALGAEGROUPS', 'SPECIES'};
tPlAbnd = grpstats(tPlAbnd(:, [GroupVars, {'Abnd'}]), GroupVars, {'sum'});
tPlAbnd.Abnd = tPlAbnd.sum_Abnd;
tPlAbnd.sum_Abnd = [];

% Unstack: species in rows, samples (DateStat) in columns
tPlAbndUns = unstack(tPlAbnd(:, {'ALGAEGROUPS', 'SPECIES', 'Abnd', 'DateStat'}), 'Abnd', 'DateStat');
tPlAbndUns = sortrows(tPlAbndUns, {'ALGAEGROUPS', 'SPECIES'});

indAbnd = 3:size(tPlAbndUns, 2);   % columns with biomass
for iCol = indAbnd
    tPlAbndUns{isnan(tPlAbndUns{:, iCol}), iCol} = 0;
end

% Keep species observed in at least two samples
ind = sum(tPlAbndUns{:, indAbnd} > 0, 2) > 1;
tPlAbndUns = tPlAbndUns(ind, :);

%% Assign environmental data to phytoplankton samples
indTempCols = {'TEMPERATURE__C', 'SALINITY__'};
indChemCols = {'pH', 'PO4_g_l', 'Ptotal_g_l', 'SiO4_g_l', 'NO2_g_l', 'NO3_g_l', ...
    'NH4_g_l', 'Ntotal_g_l', 'DIN_g_l'};

% Add environmental columns initialised with NaN
tEnvNaNs = [tTempSalDataAvgd(ones(size(tPlAbnd, 1), 1), indTempCols), ...
            tChemDataAvgd(ones(size(tPlAbnd, 1), 1), indChemCols)];
tEnvNaNs{:, :} = NaN;
tPlAbnd = [tPlAbnd, tEnvNaNs];

% For each station and sampling date take the environmental measurement
% made within +-7 days
lStations = unique(tPlAbnd.STATION);
tPlAbndEnv = tPlAbnd([], :);   % empty table with the same columns
for iS = 1:length(lStations)
    tPlAbndSel          = tPlAbnd(tPlAbnd.STATION == lStations(iS), :);
    tChemDataAvgdSel    = tChemDataAvgd(tChemDataAvgd.STATION == lStations(iS), :);
    tTempSalDataAvgdSel = tTempSalDataAvgd(tTempSalDataAvgd.STATION == lStations(iS), :);
    lDates = unique(tPlAbndSel.Date);
    for iDt = 1:size(lDates, 1)
        dtCurDate = lDates(iDt);
        indCur = tPlAbndSel.Date == dtCurDate;

        tChemDataAvgdSel.DateDiff = datenum(tChemDataAvgdSel.Date) - datenum(dtCurDate);
        ind = abs(tChemDataAvgdSel.DateDiff) < 7;
        if any(ind)
            tbl = tChemDataAvgdSel(ind, :);
            tPlAbndSel(indCur, indChemCols) = tbl(ones(sum(indCur), 1), indChemCols);
        end

        tTempSalDataAvgdSel.DateDiff = datenum(tTempSalDataAvgdSel.Date) - datenum(dtCurDate);
        ind = abs(tTempSalDataAvgdSel.DateDiff) < 7;
        if any(ind)
            tbl = tTempSalDataAvgdSel(ind, :);
            tPlAbndSel(indCur, indTempCols) = tbl(ones(sum(indCur), 1), indTempCols);
        end
    end
    tPlAbndEnv = [tPlAbndEnv; tPlAbndSel]; %#ok<AGROW>
end

%% Biomass-weighted mean environment of each species
tPlAbndEnv.DayOfYear = datenum(tPlAbndEnv.Date) - datenum(datetime(year(tPlAbndEnv.Date), 1, 1));

indEnvCols = [indTempCols, indChemCols, {'DayOfYear'}];
indEnvColsOrig = strcat(indEnvCols, '_orig');

% Keep the original values and multiply environmental values by biomass
tPlAbndEnv{:, indEnvColsOrig} = tPlAbndEnv{:, indEnvCols};
tPlAbndEnv{:, indEnvCols} = tPlAbndEnv{:, indEnvCols} .* repmat(tPlAbndEnv.Abnd, 1, numel(indEnvCols));

% Sum over all samples for each species
tPlAbndEnvSpGrpd = grpstats(tPlAbndEnv(:, [{'ALGAEGROUPS', 'SPECIES', 'Abnd'}, indEnvCols]), ...
    {'ALGAEGROUPS', 'SPECIES'}, {'sum', 'numel'});
indColsEnvMean = strcat('sum_', indEnvCols);
indColsEnvNotNaN = strcat('numel_', indEnvCols);

% grpstats returns zero if there were no measurements; replace by NaN
Env = tPlAbndEnvSpGrpd{:, indColsEnvMean};
Env(tPlAbndEnvSpGrpd{:, indColsEnvNotNaN} == 0) = NaN;
tPlAbndEnvSpGrpd{:, indColsEnvMean} = Env ./ repmat(tPlAbndEnvSpGrpd.sum_Abnd, 1, numel(indColsEnvMean));

tPlAbndUns.ID = (1:size(tPlAbndUns, 1))';
tPlAbndUnsEnv = outerjoin(tPlAbndUns, tPlAbndEnvSpGrpd, 'Keys', {'ALGAEGROUPS', 'SPECIES'}, 'Type', 'left');
tPlAbndUnsEnv = sortrows(tPlAbndUnsEnv, 'ID');

% Log-transform nutrient concentrations
indColsEnvMeanNutr2Log = {'sum_NO2_g_l', 'sum_NO3_g_l', 'sum_NH4_g_l', 'sum_Ntotal_g_l', ...
    'sum_DIN_g_l', 'sum_PO4_g_l', 'sum_Ptotal_g_l', 'sum_SiO4_g_l'};
tPlAbndUnsEnv{:, indColsEnvMeanNutr2Log} = log10(tPlAbndUnsEnv{:, indColsEnvMeanNutr2Log});

%% Optional additional plots
for iP = 1:length(PlotsID)
    switch PlotsID{iP}
        case 'VisStations'   % stations visited during the monitoring
            ShowVisStations(tPlAbnd, 10);
        otherwise
            warning('Unknown plot ID: %s', PlotsID{iP});
    end
end

%% Diffusion map
aObsCounts = sum(tPlAbndUns{:, indAbnd} > 0, 2);   % number of observations of each species
aObsMeans  = mean(tPlAbndUns{:, indAbnd}, 2);       % mean biomass of each species

% Similarity measure (see f_DM_DMit): 'Spearman' (used in the paper),
% 'Pearson', 'CosSim', 'NormzdEuc', 'StndzdEuc', 'NormzdGaus'
DM_Distance = 'Spearman';
k_max = 10;   % maximal number of nearest neighbours kept for each species
DM_DoPlots = false;   % diagnostic plots of f_DM_DMit (they use figures 1-4, which are reused below)

[ev, aEV, aS, ~, DiffDist] = f_DM_DMit(tPlAbndUns{:, indAbnd}, DM_Distance, DM_DoPlots, k_max);

% The sign of eigenvectors is arbitrary. Invert the 2nd and 3rd eigenvectors
% to reproduce the orientation obtained with MATLAB R2021 (as in the paper)
aEV(:, 2) = -aEV(:, 2);
aEV(:, 3) = -aEV(:, 3);

%% Transpose data; add date, station, functional diversity, richness and Simpson ESN
Result = DM_Calculations('Transp', -1, struct('tPlAbndUns', tPlAbndUns, 'indAbnd', indAbnd, ...
    'tStatCoords', tStatCoords, 'aEV', aEV, 'DiffDist', DiffDist));
tPlAbndTranspFreq = Result.tPlAbndTranspFreq;
indAbndTranspFreq = Result.indAbndTranspFreq;

%% Correlation between eigenvectors and environmental conditions (SI)
ColumnNames = {'sum_TEMPERATURE__C', 'sum_SALINITY__', 'sum_pH', 'sum_PO4_g_l', 'sum_Ptotal_g_l', ...
    'sum_SiO4_g_l', 'sum_NO2_g_l', 'sum_NO3_g_l', 'sum_NH4_g_l', 'sum_Ntotal_g_l', 'sum_DIN_g_l', 'sum_DayOfYear'};
ColumnNamesTitle = {'Temperature', 'salinity', 'pH', 'PO_4', 'P_{total}', 'SiO_4', ...
    'NO_2', 'NO_3', 'NH_4', 'N_{total}', 'DIN', 'Day of year'};
iFgNr = DM_Adds_Calculations('SupFig2EigV_Env_Correlation2', iFgNr, ...
    struct('aEV', aEV, 'tPlAbndUnsEnv', tPlAbndUnsEnv, 'indColsEnvMean', {indColsEnvMean}, ...
    'ColumnNames', {ColumnNames}, 'ColumnNamesTitle', {ColumnNamesTitle}));

%% Trait space coloured by environmental conditions (Fig. 3)
iFgNr = DM_Adds_Calculations('PlotEigenvectorS3DpaceWithParams', iFgNr, ...
    struct('EigVIDs', [1, 2; 1, 2; 1, 2; 1, 3], 'aEV', aEV, ...
    'EnvData', tPlAbndUnsEnv(:, {'sum_NO3_g_l', 'sum_TEMPERATURE__C', 'sum_DayOfYear', 'sum_PO4_g_l'}), ...
    'EnvDataNames', {{'log_{10}(NO_3^{-})', 'Temperature', 'Day of year', 'log_{10}(PO_4^{-})'}}));
if SaveFigures
    exportgraphics(gcf, fullfile(figDir, 'Fig3.png'), 'Resolution', 400);
end

iFgNr = DM_Adds_Calculations('PlotEigenvectorS3DpaceWithParams', iFgNr, ...
    struct('EigVIDs', [1, 2; 1, 2; 1, 3; 1, 4], 'aEV', aEV, ...
    'EnvData', tPlAbndUnsEnv(:, {'sum_SALINITY__', 'sum_NO3_g_l', 'sum_DayOfYear', 'sum_PO4_g_l'}), ...
    'EnvDataNames', {{'Salinity', 'log_{10}(NO_3^{-})', 'Day of year', 'log_{10}(PO_4^{-})'}}));

%% Biomass time series, similarity matrix and trait space (Fig. 1)
iFgNr = DM_Adds_Calculations('PlotFig1', iFgNr, struct('aEV', aEV, 'aS', aS, ...
    'SpUnStacked', tPlAbndUns{:, indAbnd}, ...
    'StationNames', {tPlAbndUns.Properties.VariableNames(indAbnd)}, 'SPECIES', {tPlAbndUns.SPECIES(:)}));
if SaveFigures
    f_exportfig_my(gcf, fullfile(figDir, 'Fig1.png'), 1);
end

iFgNr = DM_Adds_Calculations('PlotFig1Dynacom', iFgNr, struct('aEV', aEV, 'aS', aS, ...
    'SpUnStacked', tPlAbndUns{:, indAbnd}, ...
    'StationNames', {tPlAbndUns.Properties.VariableNames(indAbnd)}, 'SPECIES', {tPlAbndUns.SPECIES(:)}));

%% Functional diversity vs richness and diversity; trends in time (Fig. 4)
iFgNr = DM_Adds_Calculations('PlotFig4_FuncDivVsDiv', iFgNr, ...
    struct('tPlAbndTranspFreq', tPlAbndTranspFreq, 'indAbndTranspFreq', indAbndTranspFreq));
iFgNr = DM_Adds_Calculations('PlotFig4_FuncDivVsDivInTimeFDOnly', iFgNr, ...
    struct('tPlAbndTranspFreq', tPlAbndTranspFreq, 'indAbndTranspFreq', indAbndTranspFreq));

%% Similarity matrix, observations and mean biomass of species
subind = 1:size(aS, 1);
f_MakeFigure(1, [2274 55 1400 850]);
cm = flipud(gray);
colormap(cm);
subplot(2, 2, 1);
pcolor(subind, subind, aS(subind, subind));
colorbar
shading flat
title(['Similarity, ', DM_Distance])
f_Lbls('Species ID', 'Species ID')

subplot(2, 4, 3);
barh(subind + 0.5, aObsCounts(subind));
ylim([1, subind(end)])
f_Lbls('Observations', 'Species ID')

subplot(2, 4, 4);
barh(subind + 0.5, aObsMeans(subind));
ylim([1, subind(end)])
f_Lbls('Mean biomass', 'Species ID')

figure(2)
clf
plot(ev, '.')
title('Eigenvalues')

%% Eigenvector components and the similarity matrix sorted along them
% tSListAll lists the species with the 10 largest and 10 smallest
% components of each eigenvector (for inspection)
fg03 = f_MakeFigure(3, [100 50 1200, 1200]);
fg04 = f_MakeFigure(4, [100 50 1200, 1200]);
tSListAll = [];
ind_1st_ev = find(ev > 1e-5, 1);   % the smallest positive eigenvalue
iSubPlot = 1;
iCols = 3;
colormap(cm);
xLims = [0, size(aEV, 1)];
figure(fg03);
for iEV = ind_1st_ev:(ind_1st_ev + 5)
    if iSubPlot > 3
        iSubPlot = 1;
        figure(fg04);
    end
    [~, indEV] = sort(aEV(:, iEV), 'desc');
    subplot(3, iCols, iSubPlot);
    plot(aEV(indEV, iEV), '.')
    title(sprintf('\\lambda_{%i}=%.3g, %s', iEV, ev(iEV), DM_Distance))
    f_Lbls('Species', 'eig vect components')
    xlim(xLims);

    subplot(3, iCols, iSubPlot + iCols);
    pcolor(aS(indEV, indEV))
    shading flat
    title('Similarity matrix, sorted')
    f_Lbls('Species', 'Species')
    xlim(xLims);
    ylim(xLims);

    subplot(6, iCols, iSubPlot + 4 * iCols);
    semilogy(aObsMeans(indEV), '.g')
    f_Lbls('Species', 'Mean biom')
    xlim(xLims);

    subplot(6, iCols, iSubPlot + 5 * iCols);
    semilogy(aObsCounts(indEV), '.g')
    f_Lbls('Species', 'Observations')
    xlim(xLims);

    iSubPlot = iSubPlot + 1;

    % species with the largest components
    MaxEVComp = maxk(aEV(indEV, iEV), 10);
    indSp = aEV(:, iEV) >= MaxEVComp(end);
    tSList = tPlAbndUns(indSp, 1:2);
    tSList.Row = [];
    tSList.Lambda = ev(iEV) * ones(size(tSList.ALGAEGROUPS, 1), 1);
    tSList.Observations = aObsCounts(indSp);
    tSList.MeanBiom = aObsMeans(indSp);
    tSListAll = [tSListAll; tSList]; %#ok<AGROW>

    % species with the smallest components
    MinEVComp = mink(aEV(:, iEV), 10);
    indSp = aEV(:, iEV) <= MinEVComp(end);
    tSList = tPlAbndUns(indSp, 1:2);
    tSList.Row = [];
    tSList.Lambda = -ev(iEV) * ones(size(tSList.ALGAEGROUPS, 1), 1);
    tSList.Observations = aObsCounts(indSp);
    tSList.MeanBiom = aObsMeans(indSp);
    tSListAll = [tSListAll; tSList]; %#ok<AGROW>
end

% Sorted components of the first 16 eigenvectors
f_MakeFigure(5, [100 50 1200, 1200]);
iSubPlot = 0;
for iEV = ind_1st_ev:(ind_1st_ev + 15)
    iSubPlot = iSubPlot + 1;
    [~, indEV] = sort(aEV(:, iEV), 'desc');
    subplot(4, 4, iSubPlot);
    plot(aEV(indEV, iEV), '.')
    title(sprintf('\\lambda_{%i}=%.3g, %s', iEV, ev(iEV), DM_Distance))
    f_Lbls('Species', 'eig vect components')
end

%% Pairs of eigenvectors coloured by environmental conditions
f_MakeFigure(7, [100 50 901 992]);
atmpEnvCols = {'sum_NO3_g_l', 'sum_DayOfYear', 'sum_TEMPERATURE__C', 'sum_PO4_g_l', 'sum_pH'};
atmpEnvColsCapt = {'NO_3', 'Day of year', 'Temperature', 'PO_4', 'pH'};
aPairs = [1, 2; 1, 3; 1, 4; 1, 6];
iSubplot = 1;
for j = 1:length(atmpEnvCols)
    for i = 1:size(aPairs, 1)
        ax = subplot(length(atmpEnvCols), size(aPairs, 1), iSubplot);
        iSubplot = iSubplot + 1;
        s = scatter(aEV(:, aPairs(i, 1)), aEV(:, aPairs(i, 2)), 30, tPlAbndUnsEnv{:, atmpEnvCols{j}}, 'filled');
        s.MarkerFaceAlpha = 0.6;
        f_Lbls(sprintf('Eigenvector %i', aPairs(i, 1)), sprintf('Eigenvector %i', aPairs(i, 2)));
        caxis(quantile(tPlAbndUnsEnv{:, atmpEnvCols{j}}, [0.1, 0.8]));
    end
    pos = ax.Position;
    b = colorbar;
    b.Label.String = atmpEnvColsCapt{j};
    ax.Position = pos;
end

%% Seasonal shift of communities in the space of eigenvectors 2 and 3
f_MakeFigure(8, [100 50 1200, 1200]);

% Spring, summer and autumn samples (by day of year)
cInd = {};
cInd{1} = tPlAbndEnv.DayOfYear_orig < 100;
cInd{2} = tPlAbndEnv.DayOfYear_orig >= 100 & tPlAbndEnv.DayOfYear_orig < 250;
cInd{3} = tPlAbndEnv.DayOfYear_orig >= 250;
titles = {'Spring', 'Summer', 'Autumn'};

fPlotEllipse = @(X0, SIGMA, Clr) plot(X0(1) + SIGMA(1)*cos(-pi:0.01:pi), X0(2) + SIGMA(2)*sin(-pi:0.01:pi), ...
    'LineWidth', 2, 'Color', Clr);

colors = jet(length(cInd));
MeanXY = [];
SigmaXY = [];
for iInd = 1:length(cInd)
    subplot(2, 2, iInd);
    plot(aEV(:, 2), aEV(:, 3), '.')
    hold on
    % mean biomass of each species in this season
    tPlAbndEnvSelMean = grpstats(tPlAbndEnv(cInd{iInd}, {'ALGAEGROUPS', 'SPECIES', 'Abnd'}), ...
        {'ALGAEGROUPS', 'SPECIES'}, {'mean'});
    tTmp = outerjoin(tPlAbndUnsEnv(:, {'ALGAEGROUPS_tPlAbndUns', 'SPECIES_tPlAbndUns', 'ID'}), ...
        tPlAbndEnvSelMean(:, {'ALGAEGROUPS', 'SPECIES', 'GroupCount', 'mean_Abnd'}), ...
        'LeftKeys', {'ALGAEGROUPS_tPlAbndUns', 'SPECIES_tPlAbndUns'}, ...
        'RightKeys', {'ALGAEGROUPS', 'SPECIES'}, 'Type', 'left');
    tTmp = sortrows(tTmp, 'ID');
    scatter(aEV(:, 2), aEV(:, 3), tTmp.mean_Abnd, colors(iInd, :), 'filled')

    % biomass-weighted mean and standard deviation in trait space
    indTmp = ~isnan(tTmp.mean_Abnd);
    w = tTmp{:, {'mean_Abnd'}} ./ sum(tTmp{indTmp, {'mean_Abnd'}});
    MeanXY(iInd, :) = sum(aEV(indTmp, 2:3) .* w(indTmp, [1, 1]), 1); %#ok<SAGROW>
    SigmaXY(iInd, :) = std(aEV(indTmp, 2:3), tTmp{indTmp, {'mean_Abnd'}}); %#ok<SAGROW>

    f_Lbls('Eigenvector 2', 'Eigenvector 3');
    title(titles{iInd});
    fPlotEllipse(MeanXY(iInd, :), SigmaXY(iInd, :), colors(iInd, :)*0.8);
end

% Seasonal cycle: spring -> summer -> autumn -> spring
subplot(2, 2, 4);
plot(aEV(:, 2), aEV(:, 3), '.')
hold on
fPlotEllipse(MeanXY(1, :), SigmaXY(1, :), colors(1, :));
plot(MeanXY(1:2, 1), MeanXY(1:2, 2), 'Color', colors(1, :)*0.8, 'LineWidth', 3);
fPlotEllipse(MeanXY(2, :), SigmaXY(2, :), colors(2, :));
plot(MeanXY(2:3, 1), MeanXY(2:3, 2), 'Color', colors(2, :)*0.8, 'LineWidth', 3);
fPlotEllipse(MeanXY(3, :), SigmaXY(3, :), colors(3, :));
plot(MeanXY([3, 1], 1), MeanXY([3, 1], 2), 'Color', colors(3, :)*0.8, 'LineWidth', 3);
f_Lbls('Eigenvector 2', 'Eigenvector 3');

%% Community-weighted (CW) i-traits vs month and environment (SI)
tPlAbndTranspFreq.CWTrait1 = tPlAbndTranspFreq{:, indAbndTranspFreq} * aEV(:, 1);
% Invert the sign of CW i-trait 2 to obtain a positive correlation with temperature
tPlAbndTranspFreq.CWTrait2 = -tPlAbndTranspFreq{:, indAbndTranspFreq} * aEV(:, 2);
tPlAbndTranspFreq.CWTrait3 = tPlAbndTranspFreq{:, indAbndTranspFreq} * aEV(:, 3);

fg = figure(iFgNr);
iFgNr = iFgNr + 1;
fg.Position(3:4) = [919 713];
clf
tiledlayout('flow');

% Station categories
OffShoreStations = {'64', '65', '20', '1B'};
ExitStations = {'4'};   % exit of the Curonian Lagoon
tPlAbndTranspFreq.StationCat = tPlAbndTranspFreq.Station;
tPlAbndTranspFreq.StationCat(:) = {'on-shore'};
tPlAbndTranspFreq.StationCat(ismember(tPlAbndTranspFreq.Station, OffShoreStations)) = {'off-shore'};
tPlAbndTranspFreq.StationCat(ismember(tPlAbndTranspFreq.Station, ExitStations)) = {'exit'};

% CW i-trait 1 and mean NO3 vs month
ax = nexttile;
yyaxis left
colororder('default')
boxchart(tPlAbndTranspFreq.Month, tPlAbndTranspFreq.CWTrait1, 'GroupByColor', tPlAbndTranspFreq.StationCat)
f_Lbls('Month', 'i-trait 1, CW')
xlim([0, 13]);
ax.XTick = 1:12;
legend('Location', 'eastoutside')

tChemDataAvgd.Month = month(tChemDataAvgd.Date);
tChemDataAvgd.StationCat = tChemDataAvgd.STATION;
tChemDataAvgd.StationCat(:) = {'on-shore'};
tChemDataAvgd.StationCat(ismember(tChemDataAvgd.STATION, OffShoreStations)) = {'off-shore'};
tChemDataAvgd.StationCat(ismember(tChemDataAvgd.STATION, ExitStations)) = {'exit'};

tAvgNO3 = grpstats(tChemDataAvgd(:, {'NO3_g_l', 'Month', 'StationCat'}), {'Month', 'StationCat'}, {'@mean'});
StatCats = {'exit', 'off-shore', 'on-shore'};
yyaxis right
colororder('default')
for iSC = 1:length(StatCats)
    indSC = StatCats(iSC) == tAvgNO3.StationCat;
    p = plot(tAvgNO3.Month(indSC), tAvgNO3.Fun1_NO3_g_l(indSC), 'LineWidth', 2);
    hold on
    p.DisplayName = char(StatCats(iSC));
end
ylim([-100, 400])
ax.YTick = 0:100:400;
ylabel('Mean NO_3');

% CW i-trait 2 and mean temperature vs month
ax = nexttile;
yyaxis left
colororder('default')
boxchart(tPlAbndTranspFreq.Month, tPlAbndTranspFreq.CWTrait2, 'GroupByColor', tPlAbndTranspFreq.StationCat)
xlabel('Month');
ylabel([char(8211) '(i-trait 2), CW'])
xlim([0, 13]);
ax.XTick = 1:12;
legend('Location', 'eastoutside')

tTempSalDataAvgd.Month = month(tTempSalDataAvgd.Date);
tAvgTemp = grpstats(tTempSalDataAvgd(:, {'TEMPERATURE__C', 'Month'}), {'Month'}, {'mean'});
yyaxis right
p = plot(tAvgTemp.Month, tAvgTemp.mean_TEMPERATURE__C, 'LineWidth', 2);
p.DisplayName = 'T';
ylim([-25, 25])
ylabel('Mean temperature');
ax.YTick = -10:10:20;

AddLetters2Plots
if SaveFigures
    f_exportfig_my(gcf, fullfile(figDir, 'FigSI_EnvITraits.png'), 1);
end

%% Functional diversity aggregated by year and season
sStat = 'mean';   % statistic used for aggregation ('mean' or 'median')

% Annual mean community
tPlAbndTranspY = grpstats([tPlAbndTranspFreq(:, indAbndTranspFreq), tPlAbndTranspFreq(:, 'Year')], 'Year', sStat);
tmp_ind_abnd = 3:size(tPlAbndTranspY, 2);
FD = f_DiversityMetrics(tPlAbndTranspY{:, tmp_ind_abnd}, aEV, 'RaoDist', DiffDist, 2);
f_MakeFigure(15, [-1, -1, 1200, 1200]);
tiledlayout(3, 2, 'TileSpacing', 'compact')
nexttile
plot(tPlAbndTranspY.Year, FD, '-o', 'LineWidth', 2)

% Mean community for each year and season
tPlAbndTranspYS = grpstats([tPlAbndTranspFreq(:, indAbndTranspFreq), tPlAbndTranspFreq(:, {'Year', 'Season'})], ...
    {'Year', 'Season'}, sStat);
tmp_ind_abnd = 4:size(tPlAbndTranspYS, 2);
FD = f_DiversityMetrics(tPlAbndTranspYS{:, tmp_ind_abnd}, aEV, 'RaoDist', DiffDist, 2);
hold on
Seasons = unique(tPlAbndTranspYS.Season);
for iSsn = 1:length(Seasons)
    tmp_ind = strcmp(tPlAbndTranspYS.Season, Seasons{iSsn});
    plot(tPlAbndTranspYS.Year(tmp_ind), FD(tmp_ind), '-o');
end
legend([{'Annual'}; Seasons(:)], 'Location', 'southeast');
title('Year, season aggregated')

% Functional diversity of individual samples
FD = f_DiversityMetrics(tPlAbndTranspFreq{:, indAbndTranspFreq}, aEV, 'RaoDist', DiffDist, 2);
nexttile
boxplot(FD, tPlAbndTranspFreq.Year);
xtickangle(45)
title('Annually')

Seasons = {'Spring', 'Summer', 'Autumn', 'Winter'};
for iSsn = 1:length(Seasons)
    tmp_ind = strcmp(tPlAbndTranspFreq.Season, Seasons{iSsn});
    nexttile
    boxplot(FD(tmp_ind), tPlAbndTranspFreq.Year(tmp_ind));
    xtickangle(45)
    title(Seasons{iSsn})
end

%% Functional diversity by station and year
f_MakeFigure(12, [-1, -1, 1200, 1200]);
Stations = unique(tPlAbndTranspFreq.Station);
tiledlayout(5, 2, 'TileSpacing', 'compact')
for iSt = 1:length(Stations)
    tmp_ind = strcmp(tPlAbndTranspFreq.Station, Stations{iSt});
    nexttile
    boxplot(FD(tmp_ind), tPlAbndTranspFreq.Year(tmp_ind));
    xtickangle(45)
    title(Stations{iSt})
end

%% Linear trends of functional diversity at each station, coloured by season
f_MakeFigure(131, [-1, -1, 782 992]);
Stations = {'2', '3', '4', '5', '6', '7', '1B', '20', '64', '65'}';
tlo = tiledlayout(4, 3);
Seasons = {'Spring', 'Summer', 'Autumn'};
Clrs = f_Clrs_fresh(size(Stations, 1), 'Order', 'convenient');
Years = min(tPlAbndTranspFreq.Year):max(tPlAbndTranspFreq.Year);
iSpl = 1;
for iSt = 1:length(Stations)
    ax = nexttile;
    iSpl = iSpl + 1;
    if iSpl == 11
        % use the 10th tile for the legend
        for ist = 1:length(Seasons)
            scatter(-1, 1, 60, Clrs(ist, :), 'filled');
            hold on
        end
        xlim([0, 0.1])
        ax.XTickLabel = [];
        ax.YTickLabel = [];
        legend(Seasons, 'Location', 'northeast')
        tlo.XLabel.String = 'Year';
        tlo.YLabel.String = 'Functional diversity';
        ax = nexttile; %#ok<NASGU>
    end
    for ist = 1:length(Seasons)
        tmp_ind = strcmp(tPlAbndTranspFreq.Station, Stations{iSt}) & ...
            strcmp(tPlAbndTranspFreq.Season, Seasons{ist});
        s1 = scatter(tPlAbndTranspFreq.Year(tmp_ind), FD(tmp_ind), 90, Clrs(ist, :), 'filled');
        s1.MarkerFaceAlpha = 0.6;
        hold on
    end
    % linear regression over all seasons; * marks a significant trend
    tmp_ind = strcmp(tPlAbndTranspFreq.Station, Stations{iSt});
    ml = fitlm(tPlAbndTranspFreq.Year(tmp_ind), FD(tmp_ind));
    plot(Years, ml.predict(Years'), 'k', 'Linewidth', 2)
    if ml.Coefficients{2, 'pValue'} < 0.05
        sSign = '^*';
    else
        sSign = '';
    end
    title(sprintf('%s%s: f=%.2f + %.2GYr', Stations{iSt}, sSign, ...
        ml.Coefficients.Estimate(1), ml.Coefficients.Estimate(2)))
    ylim([0, 1.3])
    xlim([1990, 2020])
end

%% Map of the stations
f_MakeFigure(132, [-1, -1, 449 383]);
lats = [55.2 56.2];
lons = [20.6 21.4];
ShownStations = unique(tPlAbndTranspFreq(:, {'Station', 'Lat', 'Lon'}), 'rows');
ShownStations.Station = categorical(ShownStations.Station);
ShownStations.Size = ones(size(ShownStations.Station)) * 400;
colormap(f_Clrs_fresh(size(ShownStations, 1)));
geobubble(ShownStations, 'Lat', 'Lon', 'ColorVariable', 'Station', ...
    'Basemap', 'colorterrain', ...
    'BubbleColorList', f_Clrs_fresh(size(ShownStations, 1), 'Order', 'rainbow'));
geolimits(lats, lons);

%% Seasonal functional diversity and mean seasonal communities in trait space
Seasons = {'Spring', 'Summer', 'Autumn'};
tPlAbndTranspSsn = grpstats([tPlAbndTranspFreq(:, indAbndTranspFreq), tPlAbndTranspFreq(:, {'Season'})], ...
    {'Season'}, sStat);
f_MakeFigure(14, [-1, -1, 800, 800]);
subplot(2, 2, 1)
ind = ismember(tPlAbndTranspFreq.Season, Seasons);
boxplot(FD(ind), tPlAbndTranspFreq.Season(ind));
colormap(flipud(gray));
for iSn = 1:length(Seasons)
    subplot(2, 2, iSn + 1);
    ind = strcmp(tPlAbndTranspSsn.Season, Seasons{iSn});
    s1 = scatter(aEV(:, 1), aEV(:, 2), 40, tPlAbndTranspSsn{ind, 3:end}, 'filled', 'MarkerEdgeColor', 'y');
    s1.MarkerFaceAlpha = 0.8;
    title(Seasons{iSn})
end

%% Time series of temperature, NO3 and functional diversity at each station
fg09 = f_MakeFigure(9, [100 50 800, 800]);
fg11 = f_MakeFigure(11, [100 50 1266 500]);
tiledlayout(1, 2, 'TileSpacing', 'compact')
fg10 = f_MakeFigure(10, [100 50 400, 400]);

lStations = unique(tTempSalDataAvgd.STATION);
cClrs = f_Clrs_fresh(size(lStations, 1));   % one color per station
cTLeg = {};    % legend: temperature
cNLeg = {};    % legend: NO3
cLeg = {};     % legend: convex hull (not computed, panel left empty)
cLeg1 = {};    % legend: Rao functional diversity in figure 11
cRLeg = {};    % legend: Rao functional diversity in figure 9

tStatCoords.FuncDivTrend_Year = NaN(size(tStatCoords, 1), 1);
tStatCoords.FuncDivTrendRsq   = NaN(size(tStatCoords, 1), 1);
tStatCoords.FuncDivMean       = NaN(size(tStatCoords, 1), 1);
tStatCoords.FuncDivDetrendStd = NaN(size(tStatCoords, 1), 1);
dtFrom = datetime('01/10/1992', 'InputFormat', 'dd/MM/yyyy');

for iSt = 1:size(lStations, 1)
    indStatCoors = tStatCoords.StationNames == lStations(iSt);

    figure(fg09);
    % temperature
    ax1 = subplot(4, 1, 1);
    ind = tTempSalDataAvgd.STATION == lStations(iSt) & tTempSalDataAvgd.Date > dtFrom;
    [~, ~, ~, sLegend] = f_PlotAndFit(tTempSalDataAvgd.Date(ind), tTempSalDataAvgd.TEMPERATURE__C(ind), ...
        'lin', '-o', cClrs(iSt, :));
    cTLeg{iSt} = [char(lStations(iSt)) ', ', sLegend]; %#ok<SAGROW>

    % nitrate
    ax2 = subplot(4, 1, 2);
    ind = tChemDataAvgd.STATION == lStations(iSt) & tChemDataAvgd.Date > dtFrom;
    [~, ~, ~, sLegend] = f_PlotAndFit(tChemDataAvgd.Date(ind), tChemDataAvgd.NO3_g_l(ind), ...
        'lin', '-o', cClrs(iSt, :));
    cNLeg{iSt} = [char(lStations(iSt)) ', ', sLegend]; %#ok<SAGROW>

    % Rao functional diversity of each sample at this station
    indStat = tPlAbnd.STATION == lStations(iSt);
    vColNames = unique(tPlAbnd.DateStat(indStat)');
    Res = [];
    for iCN = 1:length(vColNames)
        PlAnd = tPlAbndUnsEnv{:, vColNames{iCN}};
        indPl = PlAnd > 0;
        C = strsplit(vColNames{iCN}, '_');
        SampleDate = datetime(C{2}, 'InputFormat', 'yyyyMMdd');
        Res.Dates(iCN) = SampleDate;
        Res.DayOfYear(iCN) = datenum(SampleDate) - datenum(datetime(year(SampleDate), 1, 1));
        Res.STATION(iCN) = lStations(iSt);
        Res.FuncDivRao(iCN) = f_DiversityMetrics(PlAnd(indPl), aEV(indPl, :), 'RaoDist', DiffDist(indPl, indPl));
    end
    ax3 = subplot(4, 1, 3);
    if ~isempty(Res) && length(Res.Dates) > 3
        ax4 = subplot(4, 1, 4);
        [~, ~, mFit, sLegend] = f_PlotAndFit(Res.Dates, Res.FuncDivRao, 'lin', '-o', cClrs(iSt, :));
        cRLeg{iSt} = [char(lStations(iSt)) ', ', sLegend]; %#ok<SAGROW>

        % trend (per year), mean and detrended standard deviation
        tStatCoords.FuncDivTrend_Year(indStatCoors) = mFit.Coefficients.Estimate(2) * 365;
        tStatCoords.FuncDivMean(indStatCoors) = mean(Res.FuncDivRao);
        tStatCoords.FuncDivDetrendStd(indStatCoors) = std(Res.FuncDivRao(:) - mFit.predict(datenum(Res.Dates(:))));
        tStatCoords.FuncDivTrendRsq(indStatCoors) = mFit.Rsquared.Ordinary;

        figure(fg11);
        nexttile(1);
        [~, ~, ~, sLegend] = f_PlotAndFit(Res.Dates, Res.FuncDivRao, 'lin', '.', cClrs(iSt, :));
        cLeg1{end + 1} = [char(lStations(iSt)) ', ', sLegend]; %#ok<SAGROW>

        figure(fg10);
        plot(Res.DayOfYear, Res.FuncDivRao, '.', 'Color', cClrs(iSt, :));
        hold on
    end
end
axes(ax1);
legend(cTLeg, 'Location', 'NorthEastOutside')
f_Lbls('Time', 'Temperature');

axes(ax2);
legend(cNLeg, 'Location', 'NorthEastOutside')
f_Lbls('Time', 'NO_3');

axes(ax3);
legend(cLeg, 'Location', 'NorthEastOutside')
f_Lbls('Time', 'Func div, conv hull');

axes(ax4);
legend(cRLeg, 'Location', 'NorthEastOutside')
f_Lbls('Time', 'Func div, Rao');

figure(fg11);
nexttile(1);
f_Lbls('Time', 'Func div, Rao');

%% Map of mean functional diversity and its trend at each station
f_MakeFigure(13, [-1, -1, 515 449]);
lats = [55.2 56.1];
lons = [20.0 21.6];
tStatCoords2Plot = tStatCoords(~isnan(tStatCoords.FuncDivTrend_Year), :);
tStatCoords2Plot.FuncDivTrend_Yearcat = categorical(round(tStatCoords2Plot.FuncDivTrend_Year, 3));
gb = geobubble(tStatCoords2Plot, 'WGS_84_Y_DecimalDegree', 'WGS_84_X_DecimalDegree', ...
    'SizeVariable', 'FuncDivMean', 'ColorVariable', 'FuncDivTrend_Yearcat', ...
    'Basemap', 'landcover');
geolimits(lats, lons)
gb.BubbleColorList = colormap(parula(size(unique(tStatCoords2Plot.FuncDivTrend_Yearcat), 1)));
gb.SizeLegendTitle = {'mean(Fdiv)'};
gb.ColorLegendTitle = {'\Delta Fdiv/year'};


%% Local functions

function t = SelectStations(t, lStations, dtStart)
% Keep the rows of table t measured at the stations lStations
% (and after dtStart, if it is not empty)
ind = ismember(t.STATION, lStations);
if ~isempty(dtStart)
    ind = ind & t.Date > dtStart;
end
t = t(ind, :);
end

function ShowVisStations(tPlAbnd, figID)
% Number of visits of each station and the dates of the visits.
% Stations visited at least 20 times are shown in red, at least 60 times in green.
f_MakeFigure(figID + 1, [514 85 900 800]);
tStatVsts = grpstats(tPlAbnd(:, {'Date', 'STATION'}), {'Date', 'STATION'});
tStatVsts = grpstats(tStatVsts(:, {'STATION'}), {'STATION'});
tStatVsts = sortrows(tStatVsts, 'GroupCount');
plot(tStatVsts.GroupCount, 'o');
f_Lbls('Station rank', 'Visits');

lStat60 = tStatVsts.STATION(tStatVsts.GroupCount >= 60);
lStat20 = tStatVsts.STATION(tStatVsts.GroupCount >= 20);

f_MakeFigure(figID, [514 85 900 800]);
plot(tPlAbnd.Date, tPlAbnd.STATION, 'b.');
hold on
ind = ismember(tPlAbnd.STATION, lStat20);
plot(tPlAbnd.Date(ind), tPlAbnd.STATION(ind), 'r.');
ind = ismember(tPlAbnd.STATION, lStat60);
plot(tPlAbnd.Date(ind), tPlAbnd.STATION(ind), 'g.');
end
