function s = GenerateFASTAwithMutations( varargin )
%% s = GenerateFASTAwithMutations
%do
% script to generate a fasta file w/some synonymous mutations
% for feeding to ART-SIM to generate synthetic data to test ideas
%
% LBC September 2018

NSynonVariants = 100 ; 
NWTchrs = 100 ; 
fasta_out_fn = 'sim_w_synon_subclonal_muts.fasta' ;

p = inputParser;
validPosInt = @(X) X==round(X) && X>=0 ; 
addParameter(p,'NSynonVariants',100,validPosInt);
addParameter(p,'NWTchrs',0,validPosInt);
addParameter(p,'fasta_out_fn',fasta_out_fn,@isstr);

parse(p,varargin{:}) ;
s = p.Results ; % create a struct to hold all options & args

PROJDIR = '~/Develop/Phix_mutagenesis/'  ;

% load synon / ns table
load( [ PROJDIR 'Data/DS_LanesInfo_Ilumina.mat'] );
s.D = D ; clear 'D' ; 

% load genome
s.genome = fastaread( [ PROJDIR 'ExternalData/genomeNC_001422.1.fasta']);


%% create a long fasta file w/many "chromosomes"
% each chr is a PhiX genome w/a single synonymous mutation
% optionally add N WT chrs as well
% 

warning('off','bioinfo:fastawrite:AppendToFile');
%% write NWTchrs WT chrs first
delete(s.fasta_out_fn) ; 
for I = 1:NWTchrs
    fastawrite( s.fasta_out_fn , ['WT_' num2str(I)] , s.genome.Sequence);
end

%% write NSynonVariants, each w/1 synon variant
for I = 1:NSynonVariants
    
    fastawrite( s.fasta_out_fn , ['WT_' num2str(I)] , s.genome.Sequence);
end

end