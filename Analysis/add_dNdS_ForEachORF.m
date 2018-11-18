%% Add dNdS for each gene + two shifted ORF frames (+1,+2)
%%
% 1. Generate the list of genes + their shifts with position start and
% position end.
% 2. For every ORF select part of the dataset w substitution info which
% corresponds to this ORF.
% 3. Add SynBool flag for each particular ORF
% 4. Run calculate dNdS

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


