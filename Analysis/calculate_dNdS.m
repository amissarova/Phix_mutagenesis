function dNdS = calculate_dNdS(T , proxy_error_rate_var_name , grid_error_rate , min_num_of_substitutions)
%% dNdS = calculate_dNdS(T , proxy_error_rate_var_name , grid_error_rate )
% function takes a dataset T and for any given ORF (real or artificial, defined by position_start
% and position_end) calculates dNdS. Function returns cell w/ a
% size 1xlength(grid_error_rate). For every point on a grid LLM is defined as proxy_error_rate is
% higher than certain point on grid_error_rate.
% A.M., 18.11.2018

if ~exist('min_num_of_substitutions' , 'var')
    min_num_of_substitutions = 0;
end

% select only substituions which are different from original Nt
T = T(T.IS_REF_ALLELE == 0 , :);

% add non-zero eps in order to avoid Inf values in case number of
% synonymous equals 0
eps = 0.1;

% find idx for proxy error rate
var_names_T = T.Properties.VariableNames;
idx_proxy_error_rate_var_name = find(strcmp(var_names_T , proxy_error_rate_var_name));
data = T.(idx_proxy_error_rate_var_name);
dNdS = NaN(1 , length(grid_error_rate));

% calculate total number of synonymous and nonsynonymous
dS_total = sum(T.SynBool); dS_total = dS_total + eps;
dN_total = sum(1 - T.SynBool); dN_total = dN_total + eps;
% calculate number of synonymous and nonsynonymous for every given point
% on a grid
for I = 1:length(grid_error_rate)
    idx = (data >= grid_error_rate(I));
    if sum(idx) > min_num_of_substitutions
        dS_temp = sum(T.SynBool(idx));
        dN_temp = sum(1 - T.SynBool(idx));
        dS_temp = (dS_temp+eps)/dS_total;
        dN_temp = (dN_temp+eps)/dN_total;
        dNdS(I) = dN_temp/dS_temp;
        continue
    else
        break
    end
end
