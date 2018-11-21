%% Pipeline: Add dNdS for each gene (two shifted ORF frames (+1,+2)) & each Illumina data we have
%%
% 0. Naviagte to Phix folder and pull libraries (functions in Matlab and Phix_mutagenesis)
% 1. Generate the list of genes (+ two shifted ORF frames (+1,+2)) with PositionStart and
% PositionEnd (for each Illumina data). Name: anno.
% 2. For every entry in 'anno' select corresponding genomic region from the
% corresponding Illumina dataset.
% 3. For every entry in 'anno' add SynBoolFlg (0 if particular substitution is nonsynonymous in given ORF; 1 - otherwise)
% 4. Calculate dNdS as array. Generate grid on a spectrum of proxy of error
% rate (z-score or raw Illumina error rate or something else). For every
% point on a grid select those substitutions that have corresponding proxy of
% error rate not lower than current point. Calculate number of synonymous
% and nonsynonymous amongst those substitutions and normalize each by
% total number of synonymous and nonsynonymous respectfully. Final dNdS is
% normalized number of nonsynonymous divided by normalized number of
% synonymous
% A.M., 19.11.2018

%%
% 0. naviagte to the right folder and load libraries
cd ~/Develop/Phix_mutagenesis/
addpath(genpath('~/Develop/matlab'));
addpath(genpath('~/Develop/Phix_mutagenesis'));
%%
% 1. Generate template anno w each ORF plus shifted ORFs 
% load an annotation .txt file with all Phix genes and their positions on
% the genome
anno = readtable('~/Develop/Phix_mutagenesis/ExternalData/annotationNC_001422.1.Ilumina.txt');
anno = anno(strcmp(anno.Type , 'CDS') , {'Gene' , 'PositionStart' , 'PositionEnd' });
anno.shiftORF = zeros(height(anno) , 1);
anno_ORF_shift_1 = anno; 
% add 1-NT shift ORF
anno_ORF_shift_1.shiftORF = ones(height(anno) , 1) ;
anno_ORF_shift_1.PositionStart = anno_ORF_shift_1.PositionStart + 1;
anno_ORF_shift_1.PositionEnd = anno_ORF_shift_1.PositionEnd + 1;
% add 2-NT shift ORF
anno_ORF_shift_2 = anno; 
anno_ORF_shift_2.shiftORF = repmat(2 , height(anno) , 1) ;
anno_ORF_shift_2.PositionStart = anno_ORF_shift_2.PositionStart + 2;
anno_ORF_shift_2.PositionEnd = anno_ORF_shift_2.PositionEnd + 2;
% concatenate all 3 together
anno = vertcat(anno, anno_ORF_shift_1 , anno_ORF_shift_2);

% 
% load current up-to-date table with all Illumina data
load('~/Develop/Phix_mutagenesis/Data/T.mat');
% eliminate sameNt right here to decrease memory consumption 
T = T(T.IS_REF_ALLELE == 0 , :);

% add boolean ID if substitution is from W or C strand
T.IS_ALLELE_W = char(T.NtSubstitute) == upper(char(T.NtSubstitute)) ;

% duplicate genome and add to the current position of length of genome and
% concatenate (to facilitate work with ORFs that coincide with the end of genome)
N_Phix = max(T.Position);
T_copy = T; T_copy.Position = T.Position + N_Phix; T = [T ; T_copy];
T = sortrows(T , {'Position' , 'NtSubstitute'});

%%
% Generate final ANNO_GENES table per each ORF (original and shifted) AND
% each filename
ANNO_GENES = table();
% create a list of unique Illumina datasets
unq_id = categories(T.filename);

%generate quantile-grid for z-score and raw Illumina error rate
K = 200;

data_AlleleFrequency = T.AlleleFrequency;
data_sub_prentnt_z_mode = T.sub_prentnt_z_mode;

AlleleFrequency_grid_error_rate = NaN(K , 1);
sub_prentnt_z_mode_grid_error_rate = NaN(K , 1);
for I = 1:K
    AlleleFrequency_grid_error_rate(I) = quantile(data_AlleleFrequency , (I-1)/K);
    sub_prentnt_z_mode_grid_error_rate(I) = quantile(data_sub_prentnt_z_mode , (I-1)/K);
end

%%
% for each filename generate it's DS w dNdS for each ORF (real and shifted)
for I1 = 1:length(unq_id)
    temp_anno = anno;
    temp_anno.filename = cell(height(anno) , 1); 
    temp_anno.AlleleFrequency_grid_error_rate = cell(height(anno) , 1);
    temp_anno.sub_prentnt_z_mode_grid_error_rate = cell(height(anno) , 1); 
    for I = 1:height(temp_anno)
        temp_anno.filename{I} = unq_id{I1};
        temp_anno.AlleleFrequency_grid_error_rate{I} = AlleleFrequency_grid_error_rate;
        temp_anno.sub_prentnt_z_mode_grid_error_rate{I} = sub_prentnt_z_mode_grid_error_rate;
    end
    temp_anno.AlleleFrequency_dNdS_W = cell(height(anno) , 1);
    temp_anno.AlleleFrequency_dNdS_C = cell(height(anno) , 1);
    temp_anno.sub_prentnt_z_mode_dNdS_W = cell(height(anno) , 1); 
    temp_anno.sub_prentnt_z_mode_dNdS_C = cell(height(anno) , 1); 
    for I = 1:height(temp_anno)
        fprintf('Processed file %d, %.2f.\n' , I1, I/height(anno));
        % 2. select ORF and correct filename from whole DS
        T_ORF = T(T.Position >= temp_anno.PositionStart(I) & T.Position <= temp_anno.PositionEnd(I) & ismember(T.filename , unq_id(I)), :);
        
        % 3. for given ORF add syn-bool for each position/substitution
        T_ORF = addSynOrNonsynBoolFlg(T_ORF);
        
        % 4. add dNdS-array
        % only fro reads from W-strand, z-score
        dNdS_W = calculate_dNdS(T_ORF(T_ORF.IS_ALLELE_W == 1,:) , 'sub_prentnt_z_mode' , temp_anno.sub_prentnt_z_mode_grid_error_rate{I} , 1);
        temp_anno.sub_prentnt_z_mode_dNdS_W{I} = dNdS_W;
        
        % only for reads from C-strand, z-score
        dNdS_C = calculate_dNdS(T_ORF(T_ORF.IS_ALLELE_W == 0,:) , 'sub_prentnt_z_mode' , temp_anno.sub_prentnt_z_mode_grid_error_rate{I} , 1);
        temp_anno.sub_prentnt_z_mode_dNdS_C{I} = dNdS_C;
        
        % only fro reads from W-strand, raw Illumina error rate
        dNdS_W = calculate_dNdS(T_ORF(T_ORF.IS_ALLELE_W == 1,:) , 'AlleleFrequency' , temp_anno.AlleleFrequency_grid_error_rate{I} , 1);
        temp_anno.AlleleFrequency_dNdS_W{I} = dNdS_W;
        
        % only for reads from C-strand, z-score, raw Illumina error rate
        dNdS_C = calculate_dNdS(T_ORF(T_ORF.IS_ALLELE_W == 0,:) , 'AlleleFrequency' , temp_anno.AlleleFrequency_grid_error_rate{I} , 1);
        temp_anno.AlleleFrequency_dNdS_C{I} = dNdS_C;
    end
    if isempty(ANNO_GENES)
        ANNO_GENES = temp_anno;
    else
        ANNO_GENES = vertcat(ANNO_GENES , temp_anno);
    end
    save('~/Develop/Phix_mutagenesis/Data/ANNO_GENES.mat' , 'ANNO_GENES');
end
save('~/Develop/Phix_mutagenesis/Data/ANNO_GENES.mat' , 'ANNO_GENES');

