%% load data
PDIR = '~/Develop/Phix_mutagenesis/';
addpath(genpath(PDIR));

ANNO = readtable([ PDIR 'Data/DS_MutationsInfo.tab'],'FileType','text');
ANNO = sortrows(ANNO,{'PositionNum' 'PositionNucleotide' 'PositionNucleotideSubstitute'},'ascend');


%% load all lanes from all experiments
% for each experiment, take the sums of all results from all lanes
% put all data into a single table, T

exps = {'180316' '180402'} ;
s = struct(); 
T = table();
for I = 1:numel(exps)
    s(I).exp = exps{I} ;
    s(I).ALL = table();
    s(I).lanes = dir([PDIR 'Data/lanes/' exps{I} '*.pileupCountAllAlleles.tab']); 
    for lI = 1:numel(s(I).lanes)
        s(I).ALL = vertcat( s(I).ALL , reshape_pileupCountAllAlleles(  [s(I).lanes(lI).folder filesep s(I).lanes(lI).name] ) ) ;
    end
    s(I).G = grpstats(s(I).ALL , {'chrom' 'pos' 'ref' 'alt' 'IS_REF_ALLELE'} ,'sum');
    s(I).G.exp = repmat( exps(I) , height(s(I).G) , 1);
    s(I).G.Properties.RowNames = strcat(s(I).G.Properties.RowNames,'__',s(I).G.exp) ; 
    T = vertcat( T , s(I).G) ; 
end
T = sortrows(T,{'chrom' 'pos' 'ref' 'alt' 'exp'},'ascend');
clear 's'  ; 


% per-pos mutation rate
T.mutation_rate = T.sum_allele_cov ./ T.sum_total_cov ; 

%% join with the annotation flag
Q = join(T(~T.IS_REF_ALLELE,:) , ANNO , 'LeftKey', {'pos' 'ref' 'alt'} , 'RightKey',  {'PositionNum' 'PositionNucleotide' 'PositionNucleotideSubstitute'});

%%
DATA(DATA>0.1) = 0.1 ; 
boxplot( DATA  , { strcat(Q.ref ,'->', Q.alt ) , Q.SynonymousTotalBool } )
ylabel('Mutation rate')
xlabel('Mutation')

%%
DATA = 100*Q.mutation_rate ; 

Q.GT1pct = DATA > 1 ; % great than 1% genotype abundance
Q.GT_pt1pct = DATA > 0.1 ; % great than 1% genotype abundance
Q.GT_pt01pct = DATA > 0.01 ; % great than 1% genotype abundance
G = grpstats( Q , {'ref' 'alt' 'SynonymousTotalBool' 'NonsenseTotalBool' 'MutationType' 'TriNucleotideContext'} , {'mean' 'sum'} ...
    , 'DataVars' ,{ 'GT1pct' 'GT_pt1pct' 'GT_pt01pct'})
str = arrayfun( @(I) [G.TriNucleotideContext{I}(1) G.ref{I} G.TriNucleotideContext{I}(2) '->' G.alt{I} ] , 1:height(G) ,'UniformOutput',false) ; 

idx = G.mean_GT_pt1pct > 0 ; 
bar( 100*G.mean_GT_pt1pct(idx) )
set(gca,'xtick',1:sum(idx))
set(gca,'xticklabel',str(idx));
ylabel('% of genomic positions where mutation is > 0.1%')
xlabel('Mutation')
