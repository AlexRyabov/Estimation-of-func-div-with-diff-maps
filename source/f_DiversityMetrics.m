function FD = f_DiversityMetrics(Abundance_in, Traits, Metric, Parameters, Dim)
%F_DIVERSITYMETRICS Functional and taxonomic diversity indexes.
%
%   FD = f_DiversityMetrics(Abundance, Traits, Metric, Parameters, Dim)
%
%   Inputs
%     Abundance   abundances; each column (Dim = 1, default) or each row
%                 (Dim = 2) is one community (sample)
%     Traits      species x traits matrix
%     Metric      diversity index:
%                 'Rao'        Rao's quadratic entropy with Euclidean trait
%                              distances (Botta-Dukat 2005)
%                 'RaoDist'    Rao's quadratic entropy; the distance matrix
%                              between species is given in Parameters
%                 'FD'         sum of all pairwise trait distances
%                 'FDDist'     sum of all pairwise distances given in Parameters
%                 'ConvHull'   volume of the convex hull in trait space
%                 'Hull90'     mean trait range of the species contributing
%                              95% of total abundance (at least 10 species)
%                 'Richness'   number of species with positive abundance
%                 'SimpsonESN' Simpson effective species number, 1/sum(p_i^2)
%     Parameters  distance matrix for 'RaoDist' and 'FDDist'
%     Dim         dimension along which communities are arranged (see above)
%
%   Output
%     FD          row vector with the index of each community
%
%   Reference
%     Botta-Dukat Z. (2005) Rao's quadratic entropy as a measure of
%     functional diversity based on multiple traits. J. Veg. Sci. 16, 533-540.

arguments
    Abundance_in
    Traits
    Metric
    Parameters = []
    Dim = 1
end

if Dim == 2
    Abundance_in = Abundance_in';
end

for ia = 1:size(Abundance_in, 2)
    Abundance = Abundance_in(:, ia);
    switch Metric
        case 'Rao'
            if numel(Abundance) > 2
                P = Abundance / sum(Abundance);
                D = squareform(pdist(Traits, 'euclidean'));
                FD(ia) = P' * D * P; %#ok<AGROW>
            else
                FD(ia) = 0; %#ok<AGROW>
            end
        case 'RaoDist'
            P = Abundance / sum(Abundance);
            FD(ia) = P' * Parameters * P; %#ok<AGROW>
        case 'FD'
            D = squareform(pdist(Traits, 'euclidean'));
            FD(ia) = sum(D(:)); %#ok<AGROW>
        case 'FDDist'
            FD(ia) = sum(Parameters(:)); %#ok<AGROW>
        case 'ConvHull'
            switch size(Traits, 2)
                case 2
                    DT = delaunayTriangulation(Traits(:, 1), Traits(:, 2));
                    [~, FD(ia)] = convexHull(DT); %#ok<AGROW>
                case 3
                    DT = delaunayTriangulation(Traits(:, 1), Traits(:, 2), Traits(:, 3));
                    [~, FD(ia)] = convexHull(DT); %#ok<AGROW>
                otherwise
                    [~, FD(ia)] = convhulln(Traits); %#ok<AGROW>
            end
        case 'Hull90'
            % sort abundances from large to small
            [~, ind] = sort(Abundance, 'desc');
            Abundance = Abundance(ind);
            TraitsSorted = Traits(ind, :);
            CumAbund = cumsum(Abundance);
            Idx90 = find(CumAbund > 0.95 * CumAbund(end), 1, 'first');
            if isempty(Idx90)
                Idx90 = length(CumAbund);
            else
                Idx90 = min(max(Idx90, 10), length(CumAbund));
            end
            % mean trait range of the dominant species
            TraitsSorted = TraitsSorted(1:Idx90, :);
            FD(ia) = mean(max(TraitsSorted, [], 1) - min(TraitsSorted, [], 1)); %#ok<AGROW>
        case 'Richness'
            FD(ia) = sum(Abundance > 0); %#ok<AGROW>
        case 'SimpsonESN'
            SumAbundance = sum(Abundance);
            if SumAbundance > 0
                Abundance = Abundance / SumAbundance;
                FD(ia) = 1 / sum(Abundance .* Abundance); %#ok<AGROW>
            else
                FD(ia) = 0; %#ok<AGROW>
            end
        otherwise
            error('Unknown diversity metric: %s', Metric)
    end
end
end
