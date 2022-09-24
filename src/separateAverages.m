function [gavgs, gstds] = separateAverages(comps, cols_per_day)
    [num_rows, num_cols, num_comps] = size(comps);
    gavgs = zeros(size(comps));
    gflucs = zeros(size(comps));
    gstds = zeros(num_rows, cols_per_day, num_comps);
    tdays = num_cols/cols_per_day;
    for ll = 1:num_comps
        for jj = 1:tdays
            gavgs(:, :, ll) = gavgs(:,:,ll) + ...
                repmat(comps(:, 1+(jj-1)*cols_per_day:jj*cols_per_day,ll), 1, tdays) / tdays;
        end
        gflucs(:, :, ll) = comps(:, :, ll) - gavgs(:, :, ll);
        gvar = zeros(num_rows, cols_per_day);
        for jj=1:tdays
            gvar = gvar + gflucs(:, 1+(jj-1)*cols_per_day:jj*cols_per_day, ll).^2;
        end    
        gstds(:, :, ll) = sqrt(gvar/(tdays-1));
    end
end
