function T = reshape_pileupCountAllAlleles( pileupCountAllAlleles_file_name )
%% given a pileupCountAllAlleles_file
% reshape to a three column table:
%   position alt cov
% turns it into a 'tidy' data format
% July 2018, LBC

%% load .tab file with allele counts
if ~exist(pileupCountAllAlleles_file_name,'file')
    error([pileupCountAllAlleles_file_name ' :  file not found']);
end
DATA = readtable(  pileupCountAllAlleles_file_name ,'FileType','text' ,'Delimiter','\t');

%      chrom pos ref	cov	A	C	G	T	*	-	+	X
% NC_001422.1	1	G	0	0	0	0	0	1
% NC_001422.1	2	A	1461	1461	0	0	0	0
% NC_001422.1	3	G	3673	0	0	3673	0	0
% NC_001422.1	4	T	70952	0	0	0	70952	0
% NC_001422.1	5	T	143838	0	0	0	143838	0
% NC_001422.1	6	T	176939	1	4	0	176934	0
% NC_001422.1	7	T	212733	0	8	0	212725	0
% NC_001422.1	8	A	231541	231511	1	12	17	0
%
%% reshape to one long tidy data table
ltrs = 'ACTG' ; 
T = table();
for I = 1:numel(ltrs)
    Q = table();
    Q.chrom = DATA.chrom ;
    Q.pos =  DATA.pos ;
    Q.ref = DATA.ref ; 
    Q.alt = repmat( {ltrs(I)} , height(Q),1);
    Q.total_cov = DATA.cov  ;
    Q.allele_cov = DATA.(ltrs(I)) ;
    T = vertcat( T , Q );
end

T = sortrows(T , {'pos' 'ref' 'alt'},'ascend');

% binary flag for which rows have the ref allele
r = arrayfun(@(I)T.ref{I},1:height(T));
a = arrayfun(@(I)T.alt{I},1:height(T));
T.IS_REF_ALLELE =  (r==a)' ;
end