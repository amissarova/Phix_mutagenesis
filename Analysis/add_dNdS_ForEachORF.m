%% Add dNdS for each gene + two shifted ORF frames (+1,+2)
%%
% 1. Generate the list of genes + their shifts with position start and
% position end.
% 2. For every ORF select part of the dataset w substitution info which
% corresponds to this ORF.
% 3. Add SynBool flag for each particular ORF
% 4. Run calculate dNdS

%% 0. naviagte to the right folder and load libraries
cd ~/Develop/Phix_mutagenesis/
addpath(genpath('~/Develop/matlab'));
addpath(genpath('~/Develop/Phix_mutagenesis'));

%% 1. Generate DS w each ORF plus shifted ORF

anno = dataset('file','~/Develop/Phix_mutagenesis/ExternalData/annotationNC_001422.1.Ilumina.tab');
anno = anno(strcmp(anno.Type , 'CDS') , {'PositionStart' , 'PositionEnd' , 'Name'});
anno.shiftORF = zeros(length(anno) , 1);

anno_ORF_shift_1 = anno; 
anno_ORF_shift_1.shiftORF = ones(length(anno) , 1) ;
anno_ORF_shift_1.PositionStart = anno_ORF_shift_1.PositionStart + 1;
anno_ORF_shift_1.PositionEnd = anno_ORF_shift_1.PositionEnd + 1;

anno_ORF_shift_2 = anno; 
anno_ORF_shift_2.shiftORF = repmat(2 , length(anno) , 1) ;
anno_ORF_shift_2.PositionStart = anno_ORF_shift_2.PositionStart + 2;
anno_ORF_shift_2.PositionEnd = anno_ORF_shift_2.PositionEnd + 2;

anno = [anno; anno_ORF_shift_1 ; anno_ORF_shift_2];

%% 
% in case we have a dataset w error rate info -- load it; if not -- run
% reshape_pileupCountAllAlleles
T = reshape_pileupCountAllAlleles('~/Develop/Phix_mutagenesis/Illumina_Data/180316_7001450_0407_ACC843ANXX__lane1/180316_7001450_0407_ACC8GKANXX__lane1_NoIndex_L001.pileupCountAllAlleles_StrandsSep.tab');
anno.z_grid_error_rate = cell(length(anno) , 1); 
for I = 1:length(anno)
    anno.z_grid_error_rate{I} = [-4:0.1:20];
end
N_Phix = max(T.Position);
T_copy = T; T_copy.Position = T.Position + N_Phix; T = [T ; T_copy];
T = sortrows(T , {'Position' , 'NtSubstitute'});
T = table2dataset(T);
AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');

anno.z_dNdS = cell(length(anno) , 1); 
for I = 1:length(anno)
    fprintf('Processed %.2f.\n' , I/length(anno));
    % select ORF from whole DS
    T_ORF = T(T.Position >= anno.PositionStart(I) & T.Position <= anno.PositionEnd(I), :);
    % for given ORF add syn-bool for each position/substitution
    T_ORF = addSynOrNonsynBoolFlg(T_ORF , AA_table);
    % for given ORF and given proxy of error rate (could be z-score or raw error rate from Illumina)
    % and given range for LLMs -- calculate dNdS
    z_grid_error_rate = anno.z_grid_error_rate{I};
    dNdS_vec = calculate_dNdS(T_ORF , 'sub_prentnt_z_mode' , z_grid_error_rate);
    anno.z_dNdS{I} = dNdS_vec;
end


