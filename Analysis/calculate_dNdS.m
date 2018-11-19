function dNdS = calculate_dNdS(T , proxy_error_rate_var_name , grid_error_rate )
%%
% function takes a dataset T and for any given ORF (real or artificial, defined by position_start
% and position_end) calculates dNdS. Function returns cell w/ a
% size 1xlength(grid_error_rate). For every point on a grid LLM is defined as proxy_error_rate is
% higher than certain point on grid_error_rate.
% A.M., 18.11.2018
T = T(T.IS_ALLELE_W == 1 & T.IS_REF_ALLELE == 0 , :);
var_names_T = get(T , 'VarNames');
idx_proxy_error_rate_var_name = find(strcmp(var_names_T , proxy_error_rate_var_name));
data = double(T(:  , idx_proxy_error_rate_var_name));
dNdS = NaN(1 , length(grid_error_rate));

eps = 0.1;

dS_total = sum(T.SynBool); dS_total = dS_total + eps;
dN_total = sum(1 - T.SynBool); dN_total = dN_total + eps;

for I = 1:length(grid_error_rate)
    idx = find(data >= grid_error_rate(I));
    if ~isempty(idx)
        dS_temp = sum(T.SynBool(idx));
        dN_temp = sum(1 - T.SynBool(idx));
        dS_temp = (dS_temp+eps)/dS_total;
        dN_temp = (dN_temp+eps)/dN_total;
        dNdS(I) = dN_temp/dS_temp;
    end
end
