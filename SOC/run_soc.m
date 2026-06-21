clc; clear;

% ——— Load image and form data matrix ———
fname = input('Enter image filename (e.g. ''p2.jpg''): ','s');
if ~isfile(fname)
    error('File "%s" not found in current folder: %s', fname, pwd);
end
img = imread(fname);
[f,h,k] = size(img);
x = reshape(double(img), f*h, k);   % n = f*h,  k = 3 for RGB

% ——— Parameters and SOC run ———
nk    = input('Number of clusters required: ');
if ~isscalar(nk) || nk<1 || nk~=round(nk)
    error('Number of clusters must be a positive integer.');
end

fac    = factorcal(x, nk, 1);
result = soc(x, nk, fac);

% ——— Global Silhouette ———
s   = silhouette(double(x), result.idx);
[~, GSS] = slht(s, result.idx, result.n, result.m, nk);
fprintf('\nGlobal Silhouette Score (GSI): %.4f\n\n', GSS);

% ——— Partition & Separation Indices ———
[PI, SI] = valid(result.dd, result.cc_norm, result.part.^2, nk);
fprintf('Partition Index    = %.4f\n', PI);
fprintf('Separation Index   = %.4f\n\n', SI);

% ——— Dunn Index ———
% Normalize x to [0,1] in each dimension (must match SOC’s normalization)
x_min = min(x,[],1);
x_max = max(x,[],1);
u = (x - x_min) ./ (x_max - x_min);

DI = adu(u, result.part.^2, result.cc_norm);
fprintf('Dunn Index         = %.4f\n\n', DI);

% ——— Display segmentation ———
labels = reshape(result.idx, f, h);
figure('Name','SOC Segmentation','NumberTitle','off');
imagesc(labels);
axis image off;
colormap(jet(nk));
title(sprintf('SOC Segmentation (GSI=%.3f)', GSS));