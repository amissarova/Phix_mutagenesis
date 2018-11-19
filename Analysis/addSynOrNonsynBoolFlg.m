function T = addSynOrNonsynBoolFlg(T , AA_table)
%%
% for any given ORF (real or artificial) function adds a column where each
% substition is assigned w syn-bool (1 if syn, 0 if not)
% A.M., 18.11.2018
T = sortrows(T , {'Position' , 'NtSubstitute'});
N_min = min(T.Position); N_max = max(T.Position);
if rem(N_max - N_min + 1 , 3) ~= 0
    error('Your ORF can not be converted to protein, something is wrong');
else
    T.SynBool = zeros(length(T) , 1);
    T.PositionInORF = T.Position - N_min + 1;
    for I = 1:length(T)
        if rem(T.PositionInORF(I) , 3) == 1
            nt = T.Nt(I); nt_substitute = T.NtSubstitute(I);
            post_nt = T.Nt(find(T.PositionInORF == T.PositionInORF(I) + 1 , 1)); 
            post_post_nt = T.Nt(find(T.PositionInORF == T.PositionInORF(I) + 2 , 1)); 
            str_original = strcat(char(nt) , char(post_nt) , char(post_post_nt));
            str_alternative = strcat(upper(char(nt_substitute)) , char(post_nt) , char(post_post_nt));
        elseif rem(T.PositionInORF(I) , 3) == 2
            nt = T.Nt(I); nt_substitute = T.NtSubstitute(I);
            pre_nt = T.Nt(find(T.PositionInORF == T.PositionInORF(I) - 1 , 1)); 
            post_nt = T.Nt(find(T.PositionInORF == T.PositionInORF(I) + 1 , 1)); 
            str_original = strcat(char(pre_nt) , char(nt) , char(post_nt));
            str_alternative = strcat(char(pre_nt) , upper(char(nt_substitute)) , char(post_nt));
        else
            nt = T.Nt(I); nt_substitute = T.NtSubstitute(I);
            pre_pre_nt = T.Nt(find(T.PositionInORF == T.PositionInORF(I) - 2 , 1)); 
            pre_nt = T.Nt(find(T.PositionInORF == T.PositionInORF(I) - 1 , 1)); 
            str_original = strcat(char(pre_pre_nt) , char(pre_nt) , char(nt));
            str_alternative = strcat(char(pre_pre_nt) , char(pre_nt) , upper(char(nt_substitute)));
        end
        aa_alternative = getAAfromCodon(str_alternative, AA_table);
        aa_original = getAAfromCodon(str_original , AA_table);
        T.SynBool(I) = strcmp(aa_alternative , aa_original);
    end
end

end