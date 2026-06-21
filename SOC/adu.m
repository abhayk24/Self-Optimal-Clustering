function DI = adu(u, part2, cc_norm)
    % DUNNINDEX  Dunn Index on normalized data
    %
    %   DI = dunnIndex(u, part2, cc_norm)
    %     u       n×D matrix of normalized data points
    %     part2   n×K hard‐membership matrix
    %     cc_norm K×D matrix of cluster centers
    %
    %  DI = min_{i<j} ||c_i−c_j|| / max_{m} diameter(X_m)
    %  where diameter(X_m)=max_{p,q∈C_m}||u^p−u^q||.
    
        [n,K] = size(part2);
    
        % 1) compute cluster members
        clusters = cell(1,K);
        for m = 1:K
            clusters{m} = u(logical(part2(:,m)), :);  % n_m × D
        end
    
        % 2) compute diameters
        diam = zeros(1,K);
        for m = 1:K
            pts = clusters{m};
            if size(pts,1)>1
                D = pdist(pts,'euclidean');
                diam(m) = max(D);
            else
                diam(m) = 0;
            end
        end
        maxDiam = max(diam);
    
        % 3) compute min inter‐center distance
        Dcent = pdist(cc_norm,'euclidean');  % 1×(K choose 2)
        minInter = min(Dcent);
    
        % 4) Dunn index
        DI = minInter / maxDiam;
    end