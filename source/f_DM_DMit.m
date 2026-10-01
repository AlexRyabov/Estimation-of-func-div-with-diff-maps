function [ev, aEV, aSClnd, indPos, DiffDist] = f_DM_DMit(tPlAbndUns, Norm_Method, DoPlots, k_max, Options)
%F_DM_DMIT Diffusion map of species from their abundances in samples.
%
%   [ev, aEV, aSClnd, indPos, DiffDist] = f_DM_DMit(X, Norm_Method, DoPlots, k_max)
%   [...] = f_DM_DMit(X, Norm_Method, DoPlots, k_max, Options)
%
%   Inputs
%     X            species x samples matrix (or table) of abundances
%     Norm_Method  similarity measure between species:
%                  'Spearman'   (Spearman correlation + 1)/2 (used in the paper)
%                  'Pearson'    (Pearson correlation + 1)/2
%                  'CosSim'     cosine similarity
%                  'NormzdEuc'  1/(Euclidean distance) of normalized abundances
%                  'StndzdEuc'  1/(Euclidean distance) of standardized abundances
%                  'NormzdGaus' Gaussian kernel of normalized abundances
%                  'NormzdJac'  1 - Jaccard distance
%                  'None'       X is already a similarity matrix
%     DoPlots      true to show diagnostic plots
%     k_max        maximal number of nearest neighbours kept for each species
%     Options      struct with field nboot: number of bootstrap samples used
%                  to average the similarity matrix (default 0, no bootstrap)
%
%   Outputs
%     ev        positive eigenvalues of the random-walk Laplacian, ascending
%               (at most 100)
%     aEV       corresponding eigenvectors divided by the eigenvalues,
%               i.e. species i-traits (species x eigenvectors); NaN for
%               species with zero total abundance
%     aSClnd    sparsified, symmetric similarity matrix
%     indPos    indexes of species with positive total abundance
%     DiffDist  diffusion distance between species (Euclidean distance in
%               the space of scaled eigenvectors)

if istable(tPlAbndUns)
    tPlAbndUns = tPlAbndUns{:, :};
end

if ~exist('Options', 'var')
    Options.nboot = 0;   % no bootstrapping
end
Options.bootstrap = Options.nboot > 0;
Options.SimilarityRuns = max(1, Options.nboot);   % at least one similarity matrix

% Use only species with positive total abundance
indPos = find(nansum(tPlAbndUns, 2) > 0);
tPlAbndUnsInit = tPlAbndUns;
tPlAbndUns = tPlAbndUns(indPos, :);

iSpNum = size(tPlAbndUns, 1);   % number of species
m = size(tPlAbndUns, 2);        % number of samples

%% Similarity matrix (averaged over bootstrap samples, if any)
for iSS = 1:Options.SimilarityRuns
    if Options.bootstrap
        indSample = randi(m, m, 1);   % samples drawn with replacement
        X = tPlAbndUns(:, indSample);
    else
        X = tPlAbndUns;
    end
    switch Norm_Method
        case 'CosSim'
            aNorm = repmat(sqrt(sum(X.^2, 2)), 1, size(X, 2));
            X = X ./ aNorm;
            aS = X * X';
        case 'NormzdEuc'
            aNorm = repmat(sqrt(sum(X.^2, 2)), 1, size(X, 2));
            X = X ./ aNorm;
            Dist = squareform(pdist(X, 'euclidean'));
            aS = 1 ./ Dist;
        case 'StndzdEuc'
            aMean = repmat(nanmean(X, 2), 1, size(X, 2));
            aStd = repmat(nanstd(X, 1, 2), 1, size(X, 2));
            X = (X - aMean) ./ aStd;
            Dist = squareform(pdist(X, 'euclidean'));
            aS = 1 ./ Dist;
        case 'NormzdGaus'
            aNorm = repmat(sqrt(sum(X.^2, 2)), 1, size(X, 2));
            X = X ./ aNorm;
            Dist = squareform(pdist(X, 'euclidean'));
            aS = exp(-Dist.^2);
        case 'NormzdJac'
            aNorm = repmat(sqrt(sum(X.^2, 2)), 1, size(X, 2));
            X = X ./ aNorm;
            Dist = squareform(pdist(X, 'Jaccard'));
            aS = 1 - Dist;
        case 'Pearson'
            aS = (corr(X') + 1) / 2;
        case 'Spearman'
            aS = (corr(X', 'Type', 'Spearman') + 1) / 2;
        case 'None'
            aS = X;
        otherwise
            error('Unknown similarity measure: %s', Norm_Method)
    end
    if iSS == 1
        aSTotal = aS;
    else
        aSTotal = aSTotal + aS;
    end
end
aS = aSTotal / Options.SimilarityRuns;

% Zero self-similarity
aS(1:size(aS, 1) + 1:end) = 0;

% Colour limits for the plots: 10% and 90% quantiles of similarity
MaxS = max(mink(aS(:), round(0.9 * length(aS(:)))));
MinS = min(maxk(aS(:), round(0.9 * length(aS(:)))));
subind = 1:size(aS, 1);

%% Sparsify: keep the k most similar species for each species
k = ceil(max(2, min(k_max, iSpNum)));
aSClnd = aS;
B = maxk(aS, k);
aSClnd(aS < repmat(B(end, :), size(aS, 1), 1)) = 0;
aSClnd = max(aSClnd, aSClnd');   % symmetrize

if DoPlots
    fg01 = f_MakeFigure(figure(), [10 55 1400 350]);
    cm = flipud(jet);
    cm(1, :) = [1, 1, 1];
    colormap(cm);
    subplot(1, 3, 1);
    pcolor(subind, subind, aS(subind, subind));
    colorbar
    shading flat
    title(['Sim, orig. ', Norm_Method])
    f_Lbls('Species ID', 'Species ID')

    subplot(1, 3, 2);
    pcolor(subind, subind, aSClnd(subind, subind));
    caxis([MinS, MaxS])
    colorbar
    shading flat
    title(['Sim, ' NS(k) ' sel. ', Norm_Method])
    f_Lbls('Species ID', 'Species ID')
end

%% Laplacian and its eigenvectors
% 'Gross' - random-walk Laplacian L = I - D^-1 S (used in the paper)
% 'Lafon' - Markov matrix P = D^-1 S (Coifman & Lafon)
DM_Laplacian = 'Gross';
switch DM_Laplacian
    case 'Gross'
        L = zeros(size(aSClnd));
        for i = 1:size(aSClnd, 2)
            aSClnd(i, i) = 0;
            L(i, :) = -aSClnd(i, :) / sum(aSClnd(i, :), 2);
        end
        for i = 1:size(aSClnd, 2)
            L(i, i) = 1;
        end
        [aEV, D] = eig(L);
        ev = diag(D);
        [ev, indSortEv] = sort(ev);
        aEV = aEV(:, indSortEv);

        % remove zero eigenvalues
        ind = ev > 1e-10;
        ev = ev(ind);
        aEV = aEV(:, ind);

        % keep the first 100 eigenvalues
        if length(ev) > 100
            ev = ev(1:100);
            aEV = aEV(:, 1:100);
        end
        % divide eigenvectors by eigenvalues
        aEV = repmat((1 ./ ev)', size(aEV, 1), 1) .* aEV;

    case 'Lafon'
        L = zeros(size(aSClnd));
        for i = 1:size(aSClnd, 2)
            L(i, :) = aSClnd(i, :) / sum(aSClnd(i, :), 2);
        end
        [aEV, D] = eig(L);
        ev = diag(D);
        [~, indSortEv] = sort(abs(ev), 'descend');
        ev = ev(indSortEv);
        aEV = aEV(:, indSortEv);

        % remove eigenvalues close to 1
        ind = ev < 1 - 1e-10;
        ev = ev(ind);
        aEV = aEV(:, ind);

        % multiply eigenvectors by eigenvalues
        aEV = repmat(ev(:)', size(aEV, 1), 1) .* aEV;
end

if DoPlots
    figure(fg01);
    subplot(1, 3, 3);
    PosEvInd = find(ev > 1e-6, 1);
    [~, indEV2] = sort(aEV(:, PosEvInd));   % sort along the first eigenvector
    pcolor(subind, subind, aSClnd(indEV2(subind), indEV2(subind)));
    caxis([MinS, MaxS])
    colorbar
    shading flat
    title(['Sim, v_' NS(PosEvInd) ' sorted ', Norm_Method])

    figure()
    subplot(3, 3, 1)
    plot(ev, '.')
    title('Eigenvalues')
    for iS = 2:9
        subplot(3, 3, iS)
        plot(aEV(:, iS - 1), '.')
        f_Lbls('species ID', ['ev_' NS(iS - 1)])
        title(['val=' NS(ev(iS - 1))])
    end

    figure()
    iSubplot = 1;
    for iFirstEV = 1:3
        for iSecondEV = (iFirstEV + 1):(iFirstEV + 3)
            subplot(3, 3, iSubplot);
            plot(aEV(:, iFirstEV), aEV(:, iSecondEV), '.')
            f_Lbls(['ev_' NS(iFirstEV)], ['ev_' NS(iSecondEV)])
            title(['val=' NS(ev(iFirstEV))])
            iSubplot = iSubplot + 1;
        end
    end

    figure()
    scatter3(aEV(:, 1), aEV(:, 2), aEV(:, 3), 20, aEV(:, 3), 'o', 'filled')
    f_Lbls3D('\lambda_1^{-1} v_1', '\lambda_2^{-1} v_2', '\lambda_3^{-1} v_3')
    title(['val=' NS(ev(1))])
end

%% Return results for all species (NaN for species with zero abundance)
aEV_fin = NaN(size(tPlAbndUnsInit, 1), size(aEV, 2));
aEV_fin(indPos, :) = aEV;
aEV = aEV_fin;

aSClnd_fin = NaN(size(tPlAbndUnsInit, 1), size(tPlAbndUnsInit, 1));
aSClnd_fin(indPos, indPos) = aSClnd;
aSClnd = aSClnd_fin;

if nargout >= 5
    DiffDist = squareform(pdist(aEV));
else
    DiffDist = [];
end
end
