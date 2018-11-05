function s = GenerateFASTAwithMutations( varargin )
%% s = GenerateFASTAwithMutations
%do
% script to generate a fasta file w/some synonymous mutations
% for feeding to ART-SIM to generate synthetic data to test ideas
%
% 
% eg: 
%  s = GenerateFASTAwithMutations( 'NWTchrs' , 0 , 'NSynonVariants',100 , 'fasta_out_fn' , '1pct.fasta' )
%  s = GenerateFASTAwithMutations( 'NWTchrs' , 900 , 'NSynonVariants',100 , 'fasta_out_fn' , '0.1pct.fasta' )
% LBC September 2018

NSynonVariants = 100 ; 
NWTchrs = 0 ; 
fasta_out_fn = 'sim_w_synon_subclonal_muts.fasta' ;

p = inputParser;
validPosInt = @(X) X==round(X) && X>=0 ; 
addParameter(p,'NSynonVariants',100,validPosInt);
addParameter(p,'NWTchrs',0,validPosInt);
addParameter(p,'fasta_out_fn',fasta_out_fn,@isstr);

parse(p,varargin{:}) ;
s = p.Results ; % create a struct to hold all options & args
s.start_time = datetime ; 

PROJDIR = '~/Develop/Phix_mutagenesis/'  ;

% load synon / ns table
load( [ PROJDIR 'Data/DS_LanesInfo_Ilumina.mat'] );
s.D = D ; clear 'D' ; 

% load genome
s.genome_filename = [ PROJDIR 'ExternalData/Illumina_WholeGenomeFasta/genome.fa'] ; 
s.genome = fastaread( s.genome_filename );

%% create a long fasta file w/many "chromosomes"
% each chr is a PhiX genome w/a single synonymous mutation
% optionally add N WT chrs as well
% 

warning('off','bioinfo:fastawrite:AppendToFile');
%% write NWTchrs WT chrs first
delete(s.fasta_out_fn) ; 
for I = 1:s.NWTchrs
    fastawrite( s.fasta_out_fn , ['WT_' num2str(I)] , s.genome.Sequence);
end

%% write NSynonVariants, each w/1 synon variant
s.synon_idx =  find(s.D.SynonymousTotalBool & s.D.IN_ORF);
s.chosen_synon_variants = randsample( s.synon_idx , s.NSynonVariants);
s.synon_headers = cell( s.NSynonVariants , 1) ; 
for I = 1:s.NSynonVariants
    this_seq = s.genome.Sequence ; 
    pos_to_mutate  = s.D.PositionNum(s.chosen_synon_variants(I)); 
    new_nt         = s.D.PositionNucleotideSubstitute{s.chosen_synon_variants(I)};
    prev_nt_D      = s.D.PositionNucleotide{s.chosen_synon_variants(I)};
    prev_nt_genome =  s.genome.Sequence(pos_to_mutate);
    this_seq(pos_to_mutate) = new_nt ; 
    this_seq_header = sprintf('Syn_%d_%d_%s_%s' ...
        , I , pos_to_mutate , prev_nt_genome  , new_nt ) ; 
    if prev_nt_D ~= prev_nt_genome
        error( [this_seq_header ' mismatch between D & genome']);
    end
    s.synon_headers{I} = this_seq_header ;
    fastawrite( s.fasta_out_fn , this_seq_header , this_seq );
end

s.finish_time = datetime ; 
save( regexprep(s.fasta_out_fn,'.fasta','.mat') ,  's');

end