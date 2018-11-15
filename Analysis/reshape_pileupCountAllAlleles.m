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

% % chrom	pos	ref	cov	covW	covC	A	C	G	T	a	c	g	t	*	-W	-C	+W	+C	X
% % phix	1	G	65259	38205	27054	2	0	38200	3	2	1	27050	1	0
% % phix	2	A	67538	39443	28095	39435	0	2	6	28094	0	0	1	0
% % phix	3	G	66018	38464	27550	1	0	38461	2	0	0	27549	1	0	T,T,T	t
% % phix	4	T	66246	38696	27550	2	4	2	38688	0	0	2	27548	0
% % phix	5	T	66849	39044	27805	0	1	0	39043	1	2	1	27801	0
% % phix	6	T	67742	39486	28256	2	1	2	39481	1	1	0	28254	0
% % phix	7	T	66958	38810	28148	2	3	0	38805	0	1	0	28147	0

%% reshape to one long tidy data table
%ltrs_names = { 'A' 'C' 'G' 'T' 'a' 'c' 'g' 't' 'del_base_multibase_del' 'delW' 'delC' 'insW' 'insC' 'unk_sub_indel'} ; % as of now, ignore indels
ltrs = { 'A' 'C' 'G' 'T' 'a' 'c' 'g' 't'} ;

T = table();
for I = 1:numel(ltrs)
    Q = table();
    Q.chrom = DATA.chrom ;
    Q.Position =  DATA.pos ;
    Q.Nt = DATA.ref ; 
    Q.NtSubstitute = repmat( ltrs(I) , height(Q),1);
    Q.total_cov = DATA.cov  ;
    Q.coverage_W  = DATA.covW ; 
    Q.coverage_C  = DATA.covC ; 
    Q.allele_cov = DATA.(ltrs{I}) ;
    T = vertcat( T , Q );
end

T = sortrows(T , {'Position' 'Nt' 'NtSubstitute'},'ascend');
T.Nt = categorical(T.Nt) ;
T.NtSubstitute = categorical(T.NtSubstitute) ; 

T.IS_REF_ALLELE = (T.Nt == T.NtSubstitute) | (char(T.Nt) == upper(char(T.NtSubstitute))) ; 

%% pre-nt & post-nt for the circular genome
seq = char(T.Nt(T.NtSubstitute=='A'))' ; 
Nltrs = numel(ltrs) ; % each Position has Ntltrs entries
T.PreNt = repmat( '-' , height(T) , 1); 
T.PreNt(1:Nltrs) = repmat( seq(end) , Nltrs,1);
T.PreNt( (Nltrs+1):end) = T.Nt(1:(end-Nltrs));
T.PreNt = categorical(cellstr(T.PreNt)) ; 

T.PostNt = repmat( '-' , height(T) , 1); 
T.PostNt(1:(end-Nltrs)) = T.Nt((Nltrs+1):end) ; 
T.PostNt( (end-Nltrs+1):end) = T.Nt(1:Nltrs) ;
T.PostNt = categorical(cellstr(T.PostNt)) ; 

%% calculate allele frequency
idx_is_crick = T.NtSubstitute == 'a' | T.NtSubstitute == 'c' | T.NtSubstitute == 't' | T.NtSubstitute == 'g';
T.AlleleFrequency = T.allele_cov ./ T.coverage_W ; 
T.AlleleFrequency(idx_is_crick) = T.allele_cov(idx_is_crick) ./ T.coverage_C(idx_is_crick) ; 

%% calculate zscores
id = [ char(T.PreNt) char(T.Nt) char(T.NtSubstitute)] ;  % unique ID to calc z-score
id = CharMat2CellArray(id)' ; 
uid = unique(id) ; 
T.zscore = NaN(height(T),1);
for id_I = 1:numel(uid)
    idx_into_T = ismember(id,uid{id_I}) & ~T.IS_REF_ALLELE ;
    T.zscore(idx_into_T) = zscore(T.AlleleFrequency(idx_into_T)) ; 
end


%% add filename
T.filename = repmat( {pileupCountAllAlleles_file_name} , height(T) , 1);

end