function [PI, SI] = valid(dd, cc_norm, part2, K)
    % VALID  Partition Index (PI) and Separation Index (SI) for hard‐clustering
    %
    %   [PI,SI] = valid(dd, cc_norm, part2, K)
    %     dd      n×K matrix of squared distances dd(j,m)=||u^j–c_m||^2
    %     cc_norm K×D matrix of cluster centers in normalized space
    %     part2   n×K hard‐membership matrix (0/1)
    %     K       number of clusters
    %
    %  PI = ∑_m [∑_j part2(j,m)*dd(j,m)] / [N_m * ∑_{k≠m} Dcm(k,m)]
    %  SI = [∑_m ∑_j part2(j,m)*dd(j,m)] / [n * min_{i≠j} Dcm(i,j)]
    
        [n, ~] = size(part2);
    
        % cluster sizes N_m
        N = sum(part2,1);            % 1×K
    
        % inter‐center squared distances
        Dcm = pdist2(cc_norm, cc_norm, 'squaredeuclidean');  % K×K
    
        % for PI denominator, zero out diagonal
        Dcm_noDiag = Dcm;
        Dcm_noDiag(1:K+1:end) = 0;
    
        % for SI (min_inter), ignore diagonal by setting it to Inf
        Dcm_forMin = Dcm;
        Dcm_forMin(1:K+1:end) = Inf;
    
        % numerator per cluster: sum_j part2(j,m)*dd(j,m)
        num_m = sum(part2 .* dd, 1);  % 1×K
    
        % denominator per cluster: N_m * sum_{k≠m} Dcm(k,m)
        den_m = N .* sum(Dcm_noDiag,1);  % 1×K
    
        % Partition Index
        PI = sum(num_m ./ den_m);
    
        % Separation Index
        total_num = sum(num_m);
        min_inter  = min(Dcm_forMin(:));
        SI = total_num / (n * min_inter);
    end