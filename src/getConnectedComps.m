function [corr_comps, num_comps, G] = getConnectedComps(tol, sep_scales, corr_vals, gtot)

    [num_rows, num_cols] = size(gtot);
    
    % Find the adjacency matrix and its graph
    admat = zeros(size(corr_vals));
    kpinds = log10(corr_vals) > tol;
    admat(kpinds) = 1;
    G = graph(admat);
    
    % Sort the scales so correlated groups are adjacent in the list
    [bins, binsizes] = conncomp(G);
    num_comps = length(binsizes);
    [~, Inds] = sort(bins);
    sep_scales = sep_scales(:, :, Inds);
    
    % Use adjacency graph to group correlated scales
    corr_comps = zeros(num_rows, num_cols, num_comps);
    indshft = 0;
    for ll = 1:num_comps
        corr_comps(:,:,ll) = sum(sep_scales(:,:,(1+indshft):(indshft+binsizes(ll))),3);
        indshft = indshft + binsizes(ll);
    end
end
