function iFgNrNew = DM_Adds_Calculations(PlotID, iFgNr, Data)
%DM_ADDS_CALCULATIONS Figures for DM_CompareSpec.
%
%   iFgNrNew = DM_Adds_Calculations(PlotID, iFgNr, Data)
%   draws the figure PlotID in figure number iFgNr (a figure with the name
%   PlotID is reused if it exists) and returns the next free figure number.
%
%   PlotID:
%     'PlotFig1'                          biomass time series of two species,
%                                         sorted similarity matrix and trait
%                                         space (Fig. 1)
%     'PlotFig1Dynacom'                   trait space with the relative
%                                         biomass of two samples
%     'PlotFig4_FuncDivVsDiv'             functional diversity vs richness
%                                         and Simpson ESN for each season
%     'PlotFig4_FuncDivVsDivInTimeFDOnly' functional diversity vs year
%     'PlotEigenvectorS3DpaceWithParams'  trait space coloured by
%                                         environmental conditions (Fig. 3)
%     'SupFig2EigV_Env_Correlation2'      strongest correlations between
%                                         eigenvectors and environment (SI)

if iFgNr > 0
    % Create the figure if it does not exist, otherwise clear it
    fg = findobj('Type', 'Figure', 'Name', PlotID);
    if isempty(fg)
        fg = figure(iFgNr);
        fg.Name = PlotID;
    else
        figure(fg);
    end
    iFgNrNew = iFgNr + 1;
    clf
else
    iFgNrNew = NaN;
end

switch PlotID

    case 'PlotFig1'
        fg.Position = [10 10 822 623];
        eEV = Data.aEV;
        aS = Data.aS;
        SpUnStacked = Data.SpUnStacked;
        StationNames = Data.StationNames;
        SPECIES = Data.SPECIES;
        Date = f_GetStationDate(StationNames);
        clf

        % Samples sorted by richness
        ind = SpUnStacked > 0.1;
        Richness = sum(ind, 2);
        [~, ind] = sort(Richness, 'desc');

        % Biomass time series of two selected species
        ax1 = subplot(4, 2, 1);
        iSelSample = 180;
        plot(Date, SpUnStacked(iSelSample, :), '.');
        f_Set_scales({'logy'});
        ylim([1, 3000])
        ax1.YTick = [1, 10, 100, 1000];
        ylabel('Biomass', 'FontSize', 11)
        f_FramePlotMy
        t1 = title(['\it ' SPECIES{iSelSample}]);
        t1.Position = [2767.38216129802 476.571399261784 0];
        xlim([datetime('1993-1-1'), datetime('2015-12-31')])

        ax2 = subplot(4, 2, 3);
        iSelSample = 198;
        plot(Date, SpUnStacked(iSelSample, :), '.');
        f_Set_scales({'logy'});
        ylim([1, 3000])
        ax2.YTick = [1, 10, 100, 1000];
        ylabel('Biomass', 'FontSize', 11)
        xlabel('Year', 'FontSize', 11)
        t = title(['\it ' SPECIES{iSelSample}]);
        t.Position = [2767.38216129802 476.571399261784 0];
        f_FramePlotMy
        xlim([datetime('1993-1-1'), datetime('2015-12-31')])

        % Similarity matrix sorted along the first eigenvector
        ax3 = subplot(2, 2, 2);
        [~, indEV2] = sort(eEV(:, 1));
        pcolor(aS(indEV2, indEV2));
        colormap(ax3, flipud(gray));
        cb = colorbar;
        cb.Label.String = 'Similarity';
        shading flat
        f_Lbls('Species ID, sorted', 'Species ID, sorted')
        f_FramePlotMy
        caxis([0, round(quantile(aS(aS(:) > 0), 0.9))])

        % The entire trait space (diffusion map)
        ax4 = subplot(2, 2, 3);
        cla
        iSelSample = 28;
        cm = flipud(parula);
        MrkSize = 30;
        s2 = scatter(eEV(:, 1), eEV(:, 2), MrkSize, eEV(:, 3), 'filled');
        colormap(ax4, cm)
        s2.MarkerFaceAlpha = 0.5;
        cb = colorbar;
        f_Lbls('i-trait 1, (v_1/\lambda_1)', 'i-trait 2, (v_2/\lambda_2)')
        cb.Label.String = 'i-trait 3, (v_3/\lambda_3)';
        annotation(fg, 'doublearrow', [0.253466444814199 0.345535219079313], ...
            [0.211788211788212 0.248751248751249]);
        text(0.0505199528812763, -0.0924528301886801, '\it d_{ij}', 'Rotation', 0, ...
            'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        f_FramePlotMy
        caxis([-0.59, 0.59])

        ax1.Position(3) = ax4.Position(3);
        ax1.Position(1) = ax4.Position(1);
        ax2.Position(3) = ax4.Position(3);
        ax2.Position(1) = ax4.Position(1);
        ax2.Position(2) = ax2.Position(2) + 0.01;

        % Trait space with the relative biomass of one sample
        ax5 = subplot(2, 2, 4);
        cla
        ind2 = SpUnStacked(:, ind(iSelSample)) > 0.01;
        MrkSize = 35;
        scatter(eEV(~ind2, 1), eEV(~ind2, 2), 15, [1, 1, 1]*0.75, 'filled');
        hold on
        RelBiomass = SpUnStacked(ind2, ind(iSelSample)) / sum(SpUnStacked(ind2, ind(iSelSample)));
        scatter(eEV(ind2, 1), eEV(ind2, 2), MrkSize, RelBiomass, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        colormap(ax5, cm);
        cb = colorbar;
        set(ax5, 'ColorScale', 'log')
        caxis([1e-3, 0.1])
        f_Lbls('i-trait 1, (v_1/\lambda_1)', 'i-trait 2, (v_2/\lambda_2)')
        cb.Label.String = 'Relative biomass';
        ax5.XAxis.Limits = ax4.XAxis.Limits;
        ax5.YAxis.Limits = ax4.YAxis.Limits;
        f_FramePlotMy

        [Date, Station] = f_GetStationDate(StationNames{ind(iSelSample)});
        t = title(['Station ' Station, ', ' datestr(Date, 'dd-mm-yyyy')]);
        t.Position = [-0.135878630840345 0.825578705032287 0];
        drawnow
        ax5.Position(1) = ax5.Position(1) + 0.03;
        drawnow
        ax3.Position(1) = ax5.Position(1);
        ax3.Position(3) = ax5.Position(3);
        ax3.Position(2) = ax2.Position(2);
        ax3.Position(4) = ax1.Position(2) + ax1.Position(4) - ax2.Position(2);

        AddLetters2Plots({ax1, ax3, ax4, ax5}, {'A', 'B', 'C', 'D'}, 'HShift', -0.08, 'VShift', -0.02);

    case 'PlotFig1Dynacom'
        fg.Position(3:4) = [421.6 643.2];
        iSelSamples = [28, 31];   % samples to show (ranked by richness)

        eEV = Data.aEV;
        SpUnStacked = Data.SpUnStacked;
        StationNames = Data.StationNames;
        clf

        % Samples sorted by richness
        ind = SpUnStacked > 0.1;
        Richness = sum(ind, 2);
        [~, ind] = sort(Richness, 'desc');

        for iSmp = 1:length(iSelSamples)
            nexttile
            iSelSample = iSelSamples(iSmp);
            ind2 = SpUnStacked(:, ind(iSelSample)) > 0.01;
            MrkSize = 130;
            s = scatter(eEV(~ind2, 1), eEV(~ind2, 2), MrkSize*0.4, [1, 1, 1]*0.75, 'filled', ...
                'MarkerEdgeColor', [1, 1, 1]*0.5);
            s.MarkerFaceAlpha = 0.5;
            hold on
            RelBiomass = SpUnStacked(ind2, ind(iSelSample)) / sum(SpUnStacked(ind2, ind(iSelSample)));
            s = scatter(eEV(ind2, 1), eEV(ind2, 2), MrkSize, RelBiomass, 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            s.MarkerFaceAlpha = 0.8;
            colormap(flipud(parula));
            cb = colorbar;

            ax5 = gca;
            set(ax5, 'ColorScale', 'log')
            caxis([1e-3, 0.1])
            f_Lbls('i-trait_1', 'i-trait_2')
            cb.Label.String = 'Relative biomass';
            ax5.XAxis.Label.FontSize = 14;
            ax5.YAxis.Label.FontSize = ax5.XAxis.Label.FontSize;
            cb.Label.FontSize = ax5.XAxis.Label.FontSize;
            f_FramePlotMy

            [Date, Station] = f_GetStationDate(StationNames{ind(iSelSample)});
            title(['Station ' Station, ', ' datestr(Date, 'yyyy-mm-dd')]);
            drawnow
        end

    case 'PlotFig4_FuncDivVsDiv'
        % Functional diversity vs richness (linear fit) and vs Simpson ESN
        % (Monod-type fit) for each season
        fg.Position(3:4) = [1000, 400];
        clf
        tPlAbndTranspFreq = Data.tPlAbndTranspFreq;
        Seasons = {'Winter', 'Spring', 'Summer', 'Autumn'};
        clrs = f_Clrs_fresh(4, 'Order', 'convenient');
        AT = 1;    % marker transparency
        Sz = 10;   % marker size
        xRich = linspace(0, max(tPlAbndTranspFreq.Richness));
        xESN = linspace(0, max(tPlAbndTranspFreq.SimpsonESN));
        % Monod function shifted by b(3); negative shifts are penalised
        fMonod = @(b, x) b(1) * (x - b(3)) ./ (b(2) + (x - b(3))) + 1000*heavisideNum(-b(3));
        leg1 = {};
        leg2 = {};
        for iS = 1:length(Seasons)
            ind = strcmp(tPlAbndTranspFreq.Season, Seasons{iS});
            Y = tPlAbndTranspFreq.FuncDivRao(ind);

            % functional diversity vs richness
            subplot(1, 2, 1)
            X = tPlAbndTranspFreq.Richness(ind);
            s1 = scatter(X, Y, Sz, clrs(iS, :), 'filled');
            s1.MarkerFaceAlpha = AT;
            lm = fitlm(X, Y);
            hold on
            plot_ind = xRich > 0;
            p = plot(xRich(plot_ind), lm.predict(xRich(plot_ind)'), 'Color', clrs(iS, :), 'LineWidth', 2);
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
            leg1{end + 1} = [Seasons{iS} ', R^2_{adj}=' NS(lm.Rsquared.Adjusted, 2)]; %#ok<AGROW>
            f_Lbls('Richness', 'Functional diversity')

            % functional diversity vs Simpson ESN
            subplot(1, 2, 2)
            X = tPlAbndTranspFreq.SimpsonESN(ind);
            s1 = scatter(X, Y, Sz, clrs(iS, :), 'filled');
            s1.MarkerFaceAlpha = AT;
            nlm = fitnlm(X, Y, fMonod, [1, 1, 1]);
            hold on
            plot_ind = xESN > nlm.Coefficients.Estimate(3);
            p = plot(xESN(plot_ind), nlm.predict(xESN(plot_ind)'), 'Color', clrs(iS, :), 'LineWidth', 2);
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
            leg2{end + 1} = [Seasons{iS} ', R^2_{adj}=' NS(nlm.Rsquared.Adjusted, 2)]; %#ok<AGROW>
            f_Lbls('Simpson ESN', 'Functional diversity')
        end
        subplot(1, 2, 1);
        l1 = legend(leg1, 'Location', 'Best');
        l1.Position(1:2) = [0.35399431974644 0.132312502607703];
        subplot(1, 2, 2);
        legend(leg2, 'Location', 'Best')

    case 'PlotFig4_FuncDivVsDivInTimeFDOnly'
        % Functional diversity vs year: all data and each season
        tPlAbndTranspFreq = Data.tPlAbndTranspFreq;
        Seasons = {'All', 'Spring', 'Summer', 'Autumn'};
        fg.Position(3:4) = [648 823];
        clf
        tiledlayout(4, 1, 'TileSpacing', 'compact')
        for iSsn = 1:length(Seasons)
            if strcmp(Seasons{iSsn}, 'All')
                tmp_ind = true(size(tPlAbndTranspFreq.Season));
            else
                tmp_ind = strcmp(tPlAbndTranspFreq.Season, Seasons{iSsn});
            end
            nexttile
            boxplot(tPlAbndTranspFreq.FuncDivRao(tmp_ind), tPlAbndTranspFreq.Year(tmp_ind));
            xtickangle(45)
            title(Seasons{iSsn})
            ylabel('Func. diversity')
            if iSsn == length(Seasons)
                xlabel('Year');
            end
        end
        AddLetters2Plots('HShift', -0.08, 'VShift', -0.02)

    case 'PlotEigenvectorS3DpaceWithParams'
        % Pairs of eigenvectors (i-traits) coloured by environmental conditions
        fg.Position = [10 10 688 507];
        clf
        IDs = Data.EigVIDs;          % 4 x 2 matrix: eigenvectors shown in each panel
        aEV = Data.aEV;
        EnvData = Data.EnvData;      % table with 4 environmental variables
        EnvDataNames = Data.EnvDataNames;
        for iSPl = 1:4
            ax = subplot(2, 2, iSPl);
            s1 = scatter(aEV(:, IDs(iSPl, 1)), aEV(:, IDs(iSPl, 2)), 50, EnvData{:, iSPl}, 'filled');
            % colour range
            if iSPl == 1
                caxis(quantile(EnvData{:, iSPl}, [0.1, 0.7]))
            else
                caxis(quantile(EnvData{:, iSPl}, [0.1, 0.8]))
            end
            s1.MarkerFaceAlpha = 0.6;
            hold on
            colormap(ax, parula)
            colorbar;
            f_Lbls(sprintf('i-trait %i', IDs(iSPl, 1)), sprintf('i-trait %i', IDs(iSPl, 2)))
            title(EnvDataNames{iSPl})
        end
        AddLetters2Plots(fg, {'A', 'B', 'C', 'D'}, 'HShift', -0.08, 'VShift', -0.02);

    case 'SupFig2EigV_Env_Correlation2'
        % Spearman correlations between the first eigenvectors and the
        % species mean environmental conditions. The strongest correlation
        % for each of the 9 best correlated environmental variables is shown.
        fg.Position = [10 50 800, 800];
        clf
        aEV = Data.aEV;                      % eigenvectors
        tPlAbndUnsEnv = Data.tPlAbndUnsEnv;  % table with environmental conditions
        indColsEnvMean = Data.indColsEnvMean;      % columns with environmental conditions
        ColumnNames = Data.ColumnNames;            % names of all columns
        ColumnNamesTitle = Data.ColumnNamesTitle;  % titles of all columns

        nEV = 7;      % number of eigenvectors tested
        nPanels = 9;  % number of panels shown
        aCorrCoeff = NaN(size(indColsEnvMean, 2), nEV);
        LMs = cell(size(indColsEnvMean, 2), nEV);
        for iEV = 1:nEV
            v = aEV(:, iEV);
            for iEnv = 1:length(indColsEnvMean)
                env = tPlAbndUnsEnv{:, indColsEnvMean(iEnv)};
                ind = ~isnan(v) & ~isnan(env) & ~isinf(env);
                if sum(ind) < 5
                    LMs{iEnv, iEV} = [];
                    aCorrCoeff(iEnv, iEV) = 0;
                else
                    LMs{iEnv, iEV} = fitlm(v(ind), env(ind));
                    aCorrCoeff(iEnv, iEV) = corr(v(ind), env(ind), 'Type', 'Spearman');
                end
            end
        end

        tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
        [~, indCoeff] = sort(abs(aCorrCoeff(:)), 'desc');
        Included = zeros(size(indColsEnvMean));
        iPanel = 1;
        iCoeff = 1;
        while iPanel <= nPanels
            [iEnv, iEV] = ind2sub(size(aCorrCoeff), indCoeff(iCoeff));
            iCoeff = iCoeff + 1;
            if Included(iEnv)   % each environmental variable is shown only once
                continue
            end
            Included(iEnv) = 1;
            nexttile();
            iPanel = iPanel + 1;
            v = aEV(:, iEV);
            sc1 = scatter(v, tPlAbndUnsEnv{:, indColsEnvMean(iEnv)}, 50, 'filled');
            sc1.MarkerFaceAlpha = 0.4;
            hold on
            title(sprintf('r_s=%.2f', aCorrCoeff(iEnv, iEV)))
            plot(v, LMs{iEnv, iEV}.predict(v), 'LineWidth', 2);
            ylabel(ColumnNamesTitle(strcmp(indColsEnvMean(iEnv), ColumnNames)));
            xlabel(['Trait_' NS(iEV) ', (v_' NS(iEV) '/\lambda_' NS(iEV) ')']);
            yl = ylim;
            ylim([quantile(tPlAbndUnsEnv{:, indColsEnvMean(iEnv)}, 0.01), yl(2)]);
        end

    otherwise
        error('Plot %s is not defined in DM_Adds_Calculations', PlotID)
end
end


function y = heavisideNum(x)
% Heaviside step function for numeric input (0 for x<0, 0.5 for x=0,
% 1 for x>0), identical to heaviside() of the Symbolic Math Toolbox
y = double(x > 0) + 0.5 * double(x == 0);
end
