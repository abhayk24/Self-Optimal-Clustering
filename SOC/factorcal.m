function [fac] = factorcal(x, nk, iter)
    maxIter = 15;
    factor  = ones(maxIter, nk);
    GS      = nan(1, maxIter);

    for it = iter:maxIter
        fprintf('\n--- Iteration %d ---\n', it);

        % 1) Run SOC with current factors
        result = soc(x, nk, factor(it,:));

        % 2) Print the raw thresholds δₘ
        fprintf('delta = [%s]\n', sprintf('%.4f ', result.d1));

        % 3) Global Silhouette (GSI)
        s = silhouette(double(x), result.idx);
        [S, GSi] = slht(s, result.idx, result.n, result.m, nk);
        GS(it) = GSi;
        fprintf('GSI   = %.4f\n', GSi);

        % 4) Partition & Separation Indices
        [PI, SI] = valid(result.dd, result.cc_norm, result.part.^2, nk);
        fprintf('PI    = %.4f\n', PI);
        fprintf('SI    = %.4f\n', SI);

        % 5) Dunn Index
        x_min = min(x,[],1);
        x_max = max(x,[],1);
        u     = (x - x_min) ./ (x_max - x_min);
        DI    = adu(u, result.part.^2, result.cc_norm);
        fprintf('DI    = %.4f\n', DI);

        % 6) Early‐exit checks
        if min(result.m)==0
            warning('Empty cluster encountered. Stopping early.');
            break;
        end
        if any(diff(result.d1)==0)
            warning('Duplicate δ values. Stopping early.');
            break;
        end

        % 7) Build Lagrange polynomial & solve for next δ*
        polym = lagrangepoly(result.d1, S);
        polym(end) = polym(end) - 1;      % set S(δ)=1
        rts = roots(polym);

        % pick the real δ maximizing S(δ) (or vertex if none real)
        realRts = rts(imag(rts)==0);
        if isempty(realRts)
            a = polym(1); b = polym(2);
            dmax = -b/(2*a);
        else
            % evaluate S at real roots
            Svals = polyval(polym, realRts);
            [~, idx] = min(abs(Svals - 1));
            dmax = realRts(idx);
        end

        % 8) Update factor for next iteration
        factor(it+1,:) = dmax ./ result.d1;
    end

    % pick the factor corresponding to the best GSI
    [~, bestIt] = max(GS);
    fac = factor(bestIt,:);
end