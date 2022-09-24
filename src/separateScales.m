function [sep_scales, corr_vals] = separateScales(gtot, max_lvls, wave_type)

    [nrws, nclmns] = size(gtot);
    
    sep_scales = zeros(nrws, nclmns, max_lvls);
    
    % Separate the scales
    for jj = 1:nrws
        % Get approx and detail coefficients
        [Cv1,Sv1] = wavedec(gtot(jj,:), max_lvls, wave_type);
        pterms = 0;
        for ll = 1:max_lvls+1
            % Reconstruct the signal detail by detail
            Cvp = zeros(length(Cv1), 1);
            Cvp(pterms+1:pterms+Sv1(ll)) = Cv1(pterms+1:pterms+Sv1(ll));
            sep_scales(jj,:,ll) = waverec(Cvp, Sv1, wave_type);
            pterms = pterms + Sv1(ll);
        end
    end
    
    % Compute the correlation coefficients
    corr_vals = zeros(max_lvls+1);
    for ll=1:max_lvls+1
        smatmll = scaleVariance(squeeze(sep_scales(:, 1:nclmns-1, ll)));
        smatpll = scaleVariance(squeeze(sep_scales(:, 2:nclmns, ll)));
        for jj=1:max_lvls+1
            smatmjj = scaleVariance(squeeze(sep_scales(:, 1:nclmns-1, jj)));
            smatpjj = scaleVariance(squeeze(sep_scales(:, 2:nclmns, jj)));
            corr_vals(ll, jj) = abs(mean(mean(smatmll.*smatpjj + smatpll.*smatmjj)))/2;
        end
    end
end
