function [sep_scales, corr_vals] = correlation_mat_maker(gtot, mxlvl, typ)

    [nrws, nclmns] = size(gtot);
    
    sep_scales = zeros(nrws, nclmns ,mxlvl+1);
    
    % Separate the scales
    for jj = 1:nrws
        % Get approx and detail coefficients
        [Cv1,Sv1] = wavedec(gtot(jj,:),mxlvl,typ);
        pterms = 0;
        for ll = 1:mxlvl+1
            % Reconstruct the signal detail by detail
            Cvp = zeros(length(Cv1),1);
            Cvp(pterms+1:pterms+Sv1(ll)) = Cv1(pterms+1:pterms+Sv1(ll));
            sep_scales(jj,:,ll) = waverec(Cvp,Sv1,typ);
            pterms = pterms + Sv1(ll);
        end
    end
    
    % Compute the correlation coefficients
    corr_vals = zeros(mxlvl+1);
    for ll=1:mxlvl+1
        smatmll = varscale(squeeze(sep_scales(:, 1:nclmns-1, ll)));
        smatpll = varscale(squeeze(sep_scales(:, 2:nclmns, ll)));
        for jj=1:mxlvl+1
            smatmjj = varscale(squeeze(sep_scales(:, 1:nclmns-1, jj)));
            smatpjj = varscale(squeeze(sep_scales(:, 2:nclmns, jj)));
            corr_vals(ll, jj) = abs(mean(mean(smatmll.*smatpjj + smatpll.*smatmjj)))/2;
        end
    end
end