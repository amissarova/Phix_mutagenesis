%% !! Run this pipeline once -- 
% generates DS w info about each possible substitution i.e.
% context/synonymous/nonsense/to which gene belongs/etc

cd ~/Develop/Phix_mutagenesis/
addpath(genpath('~/Develop/matlab'));

%% 1. Generate skeleton DS w all possible substitutions for each position

% read fasta file with whole sequence
PhixGenome = readSequenceFromFasta('~/Develop/Phix_mutagenesis/ExternalData/Illumina_WholeGenomeFasta/genome.fa');
PhixGenome2x = strcat(PhixGenome , PhixGenome);

% make a DS w all possible combinations
nucleotides = {'A' , 'C' , 'G' , 'T'};
nucleotides_sub = {'A' , 'C' , 'G' , 'T'};    
DS = dataset();
N = length(PhixGenome);
for I = 1:length(nucleotides)
    D = dataset();

    D.Position = NaN(N , 1);
    D.Nt = cell(N , 1);
    D.NtSubstitute = cell(N , 1);
    
    for J = 1:N
        D.Position(J) = J;
        temp_nucleotide = char(PhixGenome(J));
        D.Nt{J} = temp_nucleotide;
        D.NtSubstitute{J} = nucleotides_sub{I};
    end
    if isempty(DS)
        DS = D;
    else
        DS = [DS; D];
    end
end    
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');  
 
%% 2. add info about special features based on Illumina annotation
anno = dataset('file','~/Develop/Phix_mutagenesis/ExternalData/annotationNC_001422.1.Ilumina.tab');
idx = find(~strcmp(anno.Type , 'CDS') & ~strcmp(anno.Type , 'mRNA/exon'));
anno = anno(idx , :);
DS.SpecialFeatures = cell(length(DS) , 1);        
for I = 1:length(DS)
    idx = find(DS.Position(I) == anno.PositionStart);
    if ~isempty(idx)
        DS.SpecialFeatures{I} = anno.Type{idx};    
    end
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');  

%% 3. add info regarding to which annotated genes every position belongs
anno = dataset('file','~/Develop/Phix_mutagenesis/ExternalData/annotationNC_001422.1.Ilumina.tab');
AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');
idx = find(strcmp(anno.Type , 'CDS') );
anno = anno(idx , :);
% borders of the gene
DS.GeneBounds = cell(length(DS) , 1);
% gene names
DS.Genes = cell(length(DS) , 1);
% position of this Nt in the gene
DS.PositionInGenes = cell(length(DS) , 1);
% which AA is coded at this position
DS.AAs = cell(length(DS) , 1);
% Which position in the protein this AA has
DS.PositionInProteins = cell(length(DS) , 1);
% AA substitution after Nt at this position is mutated
DS.AASubstitutes = cell(length(DS) , 1);

N = max(DS.Position);
for I = 1:length(DS)
    fprintf('Processed %.2f.\n' , I/length(DS));
    temp_position = DS.Position(I);
    temp_GeneBounds = [];
    temp_GeneAnnotations = cell(0);
    temp_AA = cell(0);
    temp_AA_alt = cell(0);
    temp_PositionInGene = [];
    temp_PositionInProtein = [];
    for J = 1:length(anno)
        if temp_position >= anno.PositionStart(J) & temp_position <= anno.PositionEnd(J)
            temp_GeneBounds = [temp_GeneBounds; anno.PositionStart(J) anno.PositionEnd(J)];
            K = length(temp_GeneAnnotations); 
            temp_GeneAnnotations{K+1} = anno.Name{J};
            temp_PositionInGene = [temp_PositionInGene; temp_position - anno.PositionStart(J) + 1];
            
            gene = '';
            for Z = anno.PositionStart(J):anno.PositionEnd(J)
                gene = strcat(gene , PhixGenome2x(Z));
            end
            [codon , codon_num] = getCodonFromGenePosition(DS.Nt{I} , temp_position - anno.PositionStart(J) + 1 , gene);
            [codon_alt] = getCodonFromGenePosition(DS.NtSubstitute{I} , temp_position - anno.PositionStart(J) + 1 , gene);
            aa = getAAfromCodon(codon , AA_table);
            aa_alt = getAAfromCodon(codon_alt , AA_table);
            temp_AA{K+1} = aa;
            temp_AA_alt{K+1} = aa_alt;
            temp_PositionInProtein = [temp_PositionInProtein ; codon_num];
            
        elseif anno.PositionEnd(J) > N & temp_position <= anno.PositionEnd(J) - N
            temp_GeneBounds = [temp_GeneBounds; anno.PositionStart(J) anno.PositionEnd(J)];
            K = length(temp_GeneAnnotations); 
            temp_GeneAnnotations{K+1} = anno.Name{J};
            temp_PositionInGene = [temp_PositionInGene; temp_position + N - anno.PositionStart(J) + 1];
            
            gene = '';
            gene = '';
            for Z = anno.PositionStart(J):anno.PositionEnd(J)
                gene = strcat(gene , PhixGenome2x(Z));
            end
            [codon , codon_num] = getCodonFromGenePosition(DS.Nt{I} , temp_position + N - anno.PositionStart(J) + 1 , gene);
            [codon_alt] = getCodonFromGenePosition(DS.NtSubstitute{I} , temp_position + N - anno.PositionStart(J) + 1 , gene);
            aa = getAAfromCodon(codon , AA_table);
            aa_alt = getAAfromCodon(codon_alt , AA_table);
            temp_AA{K+1} = aa;
            temp_AA_alt{K+1} = aa_alt;
            temp_PositionInProtein = [temp_PositionInProtein ; codon_num];
            
        end
    end
    DS.GeneBounds{I} = temp_GeneBounds;
    DS.Genes{I} = temp_GeneAnnotations;
    DS.PositionInGenes{I} = temp_PositionInGene;
    DS.AAs{I} = temp_AA;
    DS.AASubstitutes{I} = temp_AA_alt;
    DS.PositionInProteins{I} = temp_PositionInProtein;
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');    

%% 4. Add Synonymous Bools and totalSynonymousBool

AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');
N = length(PhixGenome);
DS.SynonymousBool = cell(length(DS) , 1);
for I = 1:length(DS)
    
	fprintf('Processed %.2f.\n' , I/length(DS));
    
    PositionInGene = DS.PositionInGenes{I};
    GeneBounds = DS.GeneBounds{I};
    if ~isempty(PositionInGene)
        SynonymousBool = NaN(length(PositionInGene) , 1);
        for J = 1:length(PositionInGene)
            if GeneBounds(J,2) <= N
                GeneSeq = PhixGenome(GeneBounds(J,1):GeneBounds(J,2));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.NtSubstitute{I};
                temp_SynonymousBool = MutationSynonymousOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                SynonymousBool(J) = temp_SynonymousBool;
            else
                GeneSeq = strcat( PhixGenome(GeneBounds(J,1):N) , PhixGenome(1:GeneBounds(J,2) - N));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.NtSubstitute{I};
                temp_SynonymousBool = MutationSynonymousOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                SynonymousBool(J) = temp_SynonymousBool;
            end
        end
        DS.SynonymousBool{I} = SynonymousBool;
    end
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');      
%
DS.SynonymousTotalBool = ones(length(DS) , 1);
for I = 1:length(DS)
    temp_SynonymousBool = DS.SynonymousBool{I};
    if ~isempty(temp_SynonymousBool)
        DS.SynonymousTotalBool(I) = min(temp_SynonymousBool);
    end
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');     

%% 5. Add info if it is transition or transversion
DS.MutationType = cell(length(DS) , 1);
for I = 1:length(DS)
    DS.MutationType{I} = getMutationType(DS.Nt{I} , DS.NtSubstitute{I});
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');    

%% 6. add PreNt and PostNt
idx_first = find(DS.Position == nanmin(DS.Position) ); Nt_first = DS.Nt{idx_first(1)};
idx_last = find(DS.Position == nanmax(DS.Position) ); Nt_last = DS.Nt{idx_last(1)};

DS.PreNt = cell(length(DS) , 1);
DS.PostNt = cell(length(DS) , 1);
for I = 1:length(idx_first)
    DS.PreNt{idx_first(I)} = Nt_last;
end
for I = 1:length(idx_last)
    DS.PostNt{idx_last(I)} = Nt_first;
end
for I = 1:length(DS)
    if isempty(DS.PreNt{I})
        temp_position = DS.Position(I)-1;
        idx = find(DS.Position == temp_position , 1);
        DS.PreNt{I} = DS.Nt{idx};
    end
	if isempty(DS.PostNt{I})
        temp_position = DS.Position(I)+1;
        idx = find(DS.Position == temp_position , 1);
        DS.PostNt{I} = DS.Nt{idx};
    end
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');   
        
%% 7. Add nonsense boools and totalNonsenseBools
% add flag if mutation nonsense
AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');
N = length(PhixGenome);
DS.NonsenseBool = cell(length(DS) , 1);
for I = 1:length(DS)
    
	fprintf('Processed %.2f.\n' , I/length(DS));
    
    PositionInGene = DS.PositionInGenes{I};
    GeneBounds = DS.GeneBounds{I};
    if ~isempty(PositionInGene)
        NonsenseBool = NaN(length(PositionInGene) , 1);
        for J = 1:length(PositionInGene)
            if GeneBounds(J,2) <= N
                GeneSeq = PhixGenome(GeneBounds(J,1):GeneBounds(J,2));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.NtSubstitute{I};
                temp_NonsenseBool = MutationNonsenseOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                NonsenseBool(J) = temp_NonsenseBool;
            else
                GeneSeq = strcat( PhixGenome(GeneBounds(J,1):N) , PhixGenome(1:GeneBounds(J,2) - N));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.NtSubstitute{I};
                temp_NonsenseBool = MutationNonsenseOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                NonsenseBool(J) = temp_NonsenseBool;
            end
        end
        DS.NonsenseBool{I} = NonsenseBool;
    end
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');       
DS.NonsenseTotalBool = zeros(length(DS) , 1);
for I = 1:length(DS)
    temp_NonsenseBool = DS.NonsenseBool{I};
    if ~isempty(temp_NonsenseBool)
        DS.NonsenseTotalBool(I) = max(temp_NonsenseBool);
    end
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');     

%% 8. add InORF info and if Nt~=NtSub
DS.InORF = ones(length(DS) , 1);
DS.sameNt = zeros(length(DS) , 1);
for I = 1:length(DS)
    genes = DS.Genes{I};
    if isempty(genes)
        DS.InORF(I) = 0;
    end
    if strcmp(DS.Nt{I} , DS.NtSubstitute{I})
        DS.sameNt(I) = 1;
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS');

%% 9. Add ID of the group PreNt-Nt-NtSub
DS.PreNtNtNtSub = cell(length(DS) , 1);
for I = 1:length(DS)
    if strcmp(DS.Nt{I} , DS.NtSubstitute{I})
        DS.PreNtNtNtSub{I} = 'same';
    else
        DS.PreNtNtNtSub{I} = strcat('(' , DS.PreNt{I} , ')' , DS.Nt{I} , '->' , DS.NtSubstitute{I});
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS'); 

%% 10. add info about belonging to each gene separately
DS.A = zeros(length(DS) , 1); DS.A2 = zeros(length(DS) , 1); DS.B = zeros(length(DS) , 1);
DS.C = zeros(length(DS) , 1); DS.D = zeros(length(DS) , 1); DS.E = zeros(length(DS) , 1);
DS.F = zeros(length(DS) , 1); DS.G = zeros(length(DS) , 1); DS.H = zeros(length(DS) , 1);
DS.J = zeros(length(DS) , 1); DS.K = zeros(length(DS) , 1);

DS.A_Pos = NaN(length(DS) , 1); DS.A2_Pos = NaN(length(DS) , 1); DS.B_Pos = NaN(length(DS) , 1);
DS.C_Pos = NaN(length(DS) , 1); DS.D_Pos = NaN(length(DS) , 1); DS.E_Pos = NaN(length(DS) , 1);
DS.F_Pos = NaN(length(DS) , 1); DS.G_Pos = NaN(length(DS) , 1); DS.H_Pos = NaN(length(DS) , 1);
DS.J_Pos = NaN(length(DS) , 1); DS.K_Pos = NaN(length(DS) , 1);

var_names = get(DS , 'VarNames');

for I = 1:length(DS)
    genes = DS.Genes{I};
    positions_in_genes = DS.PositionInGenes{I};
    if ~isempty(genes)
        for J = 1:length(genes)
            idx = find(strcmp(genes{J} , var_names));
            DS{I , idx} = 1;
            idx = find(strcmp(strcat(genes{J}, '_Pos') , var_names));
            DS{I , idx} = positions_in_genes(J);
        end
    end
end
DS.SynBoolPerPosInAA = cell(length(DS) , 1);
AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');
N = nanmax(DS.Position);
for I = 1:length(DS)
    SynBoolPerPosinAA = [0 0 0];
    temp_position = DS.Position(I);
    if temp_position == 1
        idx1 = find(DS.Position == N-1 , 1); 
        idx2 = find(DS.Position == N , 1); 
        idx3 = find(DS.Position == 1 , 1); 
        idx4 = find(DS.Position == 2 , 1); 
        idx5 = find(DS.Position == 3 , 1);
    elseif temp_position == 2
        idx1 = find(DS.Position == N , 1); 
        idx2 = find(DS.Position == 1 , 1); 
        idx3 = find(DS.Position == 2 , 1); 
        idx4 = find(DS.Position == 3 , 1); 
        idx5 = find(DS.Position == 4 , 1);
    elseif temp_position == N-1
        idx1 = find(DS.Position == N-3 , 1); 
        idx2 = find(DS.Position == N-2 , 1); 
        idx3 = find(DS.Position == N-1 , 1); 
        idx4 = find(DS.Position == N , 1); 
        idx5 = find(DS.Position == 1 , 1);
    elseif temp_position == N
        idx1 = find(DS.Position == N-2 , 1); 
        idx2 = find(DS.Position == N-1 , 1); 
        idx3 = find(DS.Position == N , 1); 
        idx4 = find(DS.Position == 1 , 1); 
        idx5 = find(DS.Position == 2 , 1);
    else
        idx1 = find(DS.Position == temp_position-2 , 1); 
        idx2 = find(DS.Position == temp_position-1 , 1); 
        idx3 = find(DS.Position == temp_position , 1); 
        idx4 = find(DS.Position == temp_position+1 , 1); 
        idx5 = find(DS.Position == temp_position+2 , 1);
    end
    
    str3 = strcat(DS.Nt{idx1} , DS.Nt{idx2} , DS.Nt{idx3});
	str3_alt = strcat(DS.Nt{idx1} , DS.Nt{idx2} , DS.NtSubstitute{idx3});
    AA3 = getAAfromCodon(str3 , AA_table); AA3_alt = getAAfromCodon(str3_alt , AA_table); 
    if strcmp(AA3 , AA3_alt) == 1
        SynBoolPerPosinAA(3) = 1;
    end
    
    str2 = strcat(DS.Nt{idx2} , DS.Nt{idx3} , DS.Nt{idx4});
	str2_alt = strcat(DS.Nt{idx2} , DS.NtSubstitute{idx3} , DS.Nt{idx4});
    AA2 = getAAfromCodon(str2 , AA_table); AA2_alt = getAAfromCodon(str2_alt , AA_table); 
    if strcmp(AA2 , AA2_alt) == 1
        SynBoolPerPosinAA(2) = 1;
    end
    
    str1 = strcat(DS.NtSubstitute{idx3} , DS.Nt{idx4} , DS.Nt{idx5});
	str1_alt = strcat(DS.Nt{idx3} , DS.Nt{idx4} , DS.Nt{idx5});
    AA1 = getAAfromCodon(str1 , AA_table); AA1_alt = getAAfromCodon(str1_alt , AA_table); 
    if strcmp(AA1 , AA1_alt) == 1
        SynBoolPerPosinAA(1) = 1;
    end
    
    DS.SynBoolPerPosInAA{I} = SynBoolPerPosinAA;   
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat' , 'DS'); 

%% cosmetic changes and convert for the table
load('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat');
D = DS;
D.NtSubstitute = lower(D.NtSubstitute);
DS = vertcat(DS , D);
DS = sortrows(DS , {'Position' , 'NtSubstitute'});
DS.PreNtNtNtSub = categorical(DS.PreNtNtNtSub);
DS.Nt = categorical(DS.Nt);
DS.PostNt = categorical(DS.PostNt);
DS.PreNt = categorical(DS.PreNt);
DS.NtSubstitute = categorical(DS.NtSubstitute);
DS = DS(: , [1:42]);
%%
DS = dataset2table(DS);
ANNO_SUBSTITUTION = DS;
save('~/Develop/Phix_Mutagenesis/Data/ANNO_SUBSTITUTION.mat' , 'ANNO_SUBSTITUTION');











