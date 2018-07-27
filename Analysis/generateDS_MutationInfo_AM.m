cd ~/Develop/Phix_mutagenesis/
addpath(genpath('~/Develop/matlab'));

%% Pipeline to generate DS w/ all possible mutations and their annotations
% read fasta file with whole sequence
PhixGenome = readSequenceFromFasta('~/Develop/Phix_mutagenesis/ExternalData/genomeNC_001422.1.fasta');
PhixGenome2x = strcat(PhixGenome , PhixGenome);

% make a DS w all possible combinations
nucleotides_cell = {'A' , 'C' , 'G' , 'T'};
nucleotides_sub_cell = cell( length(nucleotides_cell) , length(nucleotides_cell) - 1);
for I = 1:length(nucleotides_cell)
    rest_nucleotides = setdiff(nucleotides_cell , nucleotides_cell{I});
    nucleotides_sub_cell(I,:) = rest_nucleotides;
end    
DS = dataset();
N = length(PhixGenome);
for I = 1:length(nucleotides_cell) - 1
    D = dataset();

    D.PositionNum = NaN(N , 1);
    D.PositionNucleotide = cell(N , 1);
    D.PositionNucleotideSubstitute = cell(N , 1);
    
    for J = 1:N
        D.PositionNum(J) = J;
        temp_nucleotide = char(PhixGenome(J));
        D.PositionNucleotide{J} = temp_nucleotide;
        idx = find(strcmp(nucleotides_cell , temp_nucleotide));
        D.PositionNucleotideSubstitute{J} = nucleotides_sub_cell{idx, I};
    end
    if isempty(DS)
        DS = D;
    else
        DS = [DS; D];
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');      
        
% add info reagarding a nucleotide being a common SNP or ARS
anno = dataset('file','~/Develop/Phix_mutagenesis/ExternalData/annotationNC_001422.1.tab');
idx = find(~strcmp(anno.Type , 'CDS') & ~strcmp(anno.Type , 'mRNA/exon'));
anno = anno(idx , :);
DS.SpecialFeatures = cell(length(DS) , 1);        
for I = 1:length(DS)
    idx = find(DS.PositionNum(I) == anno.PositionStart);
    if ~isempty(idx)
        DS.SpecialFeatures{I} = anno.Type{idx};    
    end
end
%
% add info regarding to which annotated genes position belongs
anno = dataset('file','~/Develop/Phix_mutagenesis/ExternalData/annotationNC_001422.1.tab');
AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');
idx = find(strcmp(anno.Type , 'CDS') );
anno = anno(idx , :);
DS.GeneBounds = cell(length(DS) , 1);
DS.GeneAnnotations = cell(length(DS) , 1);
DS.PositionInGene = cell(length(DS) , 1);

DS.AA = cell(length(DS) , 1);
DS.PositionInProtein = cell(length(DS) , 1);
DS.AASubstitutes = cell(length(DS) , 1);

N = max(DS.PositionNum);
for I = 1:length(DS)
    
    temp_position = DS.PositionNum(I);
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
            [codon , codon_num] = getCodonFromGenePosition(DS.PositionNucleotide{I} , temp_position - anno.PositionStart(J) + 1 , gene);
            [codon_alt] = getCodonFromGenePosition(DS.PositionNucleotideSubstitute{I} , temp_position - anno.PositionStart(J) + 1 , gene);
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
            [codon , codon_num] = getCodonFromGenePosition(DS.PositionNucleotide{I} , temp_position + N - anno.PositionStart(J) + 1 , gene);
            [codon_alt] = getCodonFromGenePosition(DS.PositionNucleotideSubstitute{I} , temp_position + N - anno.PositionStart(J) + 1 , gene);
            aa = getAAfromCodon(codon , AA_table);
            aa_alt = getAAfromCodon(codon_alt , AA_table);
            temp_AA{K+1} = aa;
            temp_AA_alt{K+1} = aa_alt;
            temp_PositionInProtein = [temp_PositionInProtein ; codon_num];
            
        end
    end
    DS.GeneBounds{I} = temp_GeneBounds;
    DS.GeneAnnotations{I} = temp_GeneAnnotations;
    DS.PositionInGene{I} = temp_PositionInGene;
    DS.AA{I} = temp_AA;
    DS.AASubstitutes{I} = temp_AA_alt;
    DS.PositionInProtein{I} = temp_PositionInProtein;
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');    

%
AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');
N = length(PhixGenome);
DS.SynonymousBool = cell(length(DS) , 1);
for I = 1:length(DS)
    
	fprintf('Processed %.2f.\n' , I/length(DS));
    
    PositionInGene = DS.PositionInGene{I};
    GeneBounds = DS.GeneBounds{I};
    if ~isempty(PositionInGene)
        SynonymousBool = NaN(length(PositionInGene) , 1);
        for J = 1:length(PositionInGene)
            if GeneBounds(J,2) <= N
                GeneSeq = PhixGenome(GeneBounds(J,1):GeneBounds(J,2));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.PositionNucleotideSubstitute{I};
                temp_SynonymousBool = MutationSynonymousOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                SynonymousBool(J) = temp_SynonymousBool;
            else
                GeneSeq = strcat( PhixGenome(GeneBounds(J,1):N) , PhixGenome(1:GeneBounds(J,2) - N));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.PositionNucleotideSubstitute{I};
                temp_SynonymousBool = MutationSynonymousOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                SynonymousBool(J) = temp_SynonymousBool;
            end
        end
        DS.SynonymousBool{I} = SynonymousBool;
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');      
%
DS.SynonymousTotalBool = ones(length(DS) , 1);
for I = 1:length(DS)
    temp_SynonymousBool = DS.SynonymousBool{I};
    if ~isempty(temp_SynonymousBool)
        DS.SynonymousTotalBool(I) = min(temp_SynonymousBool);
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');     

%
DS.MutationType = cell(length(DS) , 1);
for I = 1:length(DS)
    DS.MutationType{I} = getMutationType(DS.PositionNucleotide{I} , DS.PositionNucleotideSubstitute{I});
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');    

%
DS.TriNucleotideContext = cell(length(DS) , 1);

DS.TriNucleotideContext{1} = strcat(DS.PositionNucleotide{length(DS)} , DS.PositionNucleotide{2});
DS.TriNucleotideContext{length(DS)} = strcat(DS.PositionNucleotide{length(DS)-1} , DS.PositionNucleotide{1});

for I = 2:length(DS)-1
	DS.TriNucleotideContext{I} = strcat(DS.PositionNucleotide{I-1} , DS.PositionNucleotide{I+1});
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');   
        
% add flag if mutation nonsense
AA_table = dataset('file' , '~/Develop/Phix_mutagenesis/ExternalData/AA_table.tab');
N = length(PhixGenome);
DS.NonsenseBool = cell(length(DS) , 1);
for I = 1:length(DS)
    
	%fprintf('Processed %.2f.\n' , I/length(DS));
    
    PositionInGene = DS.PositionInGene{I};
    GeneBounds = DS.GeneBounds{I};
    if ~isempty(PositionInGene)
        NonsenseBool = NaN(length(PositionInGene) , 1);
        for J = 1:length(PositionInGene)
            if GeneBounds(J,2) <= N
                GeneSeq = PhixGenome(GeneBounds(J,1):GeneBounds(J,2));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.PositionNucleotideSubstitute{I};
                temp_NonsenseBool = MutationNonsenseOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                NonsenseBool(J) = temp_NonsenseBool;
            else
                GeneSeq = strcat( PhixGenome(GeneBounds(J,1):N) , PhixGenome(1:GeneBounds(J,2) - N));
                PositionNum = PositionInGene(J);
                PositionNucleotideSubstitute = DS.PositionNucleotideSubstitute{I};
                temp_NonsenseBool = MutationNonsenseOrNot(GeneSeq , PositionNum, PositionNucleotideSubstitute , AA_table);
                NonsenseBool(J) = temp_NonsenseBool;
            end
        end
        DS.NonsenseBool{I} = NonsenseBool;
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');     

%
DS.NonsenseTotalBool = zeros(length(DS) , 1);
for I = 1:length(DS)
    temp_NonsenseBool = DS.NonsenseBool{I};
    if ~isempty(temp_NonsenseBool)
        DS.NonsenseTotalBool(I) = max(temp_NonsenseBool);
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS');   

%% add a score for conservation
load('~/Develop/Phix_Mutagenesis/Data/DS_PhylogeneticInfo.mat' );   
% add vector of alternatives

DS.ConservationScore = cell(length(DS) , 1);
for I = 1:length(DS)
    
    if DS.SynonymousTotalBool(I) == 1
        DS.ConservationScore{I} = 'S';
    else
        AA = DS.AA{I};
        AASubstitutes = DS.AASubstitutes{I};
        PositionInProtein = DS.PositionInProtein{I};
        Genes = DS.GeneAnnotations{I};
        K = length(AA);
        count = 0;
        for J = 1:K
            if strcmp(AA{J} , AASubstitutes{J})
                count = count + 1;
                continue
            else 
                idx = find(strcmp(D.Gene , Genes{J}) & D.AAPosition == PositionInProtein(J));
                if isempty(idx) & strcmp(AA{J} , 'STOP') & ~strcmp(AASubstitutes{J} , 'STOP')
                    DS.ConservationScore{I} = 'N';
                    break
                elseif ~isempty(find(strcmp(D.AASpectrum{idx} , '-')))
                    DS.ConservationScore{I} = 'X';
                    break
                elseif isempty(find(strcmp(D.AASpectrum{idx} , AASubstitutes{J})))
                    DS.ConservationScore{I} = 'N';
                    break
                else
                    count = count + 1;
                    continue
                end
            end
            
        end
        if count == K
            DS.ConservationScore{I} = 'C';
        end
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.mat' , 'DS'); 


%%
load('~/Develop/Phix_mutagenesis/Data/DS_MutationsInfo.mat');
D = DS(: , {'PositionNum', 'PositionNucleotide' , 'PositionNucleotideSubstitute',...
    'SynonymousTotalBool' , 'MutationType' , 'TriNucleotideContext',...
    'NonsenseTotalBool'});
export(D , 'file' , '~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo.tab');







        