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
    Q.coverage_total = DATA.cov  ;
    Q.coverage_W  = DATA.covW ; 
    Q.coverage_C  = DATA.covC ; 
    Q.AlleleCount = DATA.(ltrs{I}) ;
    T = vertcat( T , Q );
end

T = sortrows(T , {'Position' 'Nt' 'NtSubstitute'},'ascend');
T.Nt = categorical(T.Nt) ;
T.NtSubstitute = categorical(T.NtSubstitute) ; 

T.IS_REF_ALLELE = logical( (T.Nt == T.NtSubstitute) | (char(T.Nt) == upper(char(T.NtSubstitute))) ); 
T.IS_ALLELE_W = logical( char(T.NtSubstitute) == upper(char(T.NtSubstitute)) ) ; 

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
T.AlleleFrequency = T.AlleleCount ./ T.coverage_W ; 
T.AlleleFrequency(idx_is_crick) = T.AlleleCount(idx_is_crick) ./ T.coverage_C(idx_is_crick) ; 

%% Categorical values for substitution IDs
T.sub_prentnt = categorical( CharMat2CellArray([ char(T.NtSubstitute) char(T.PreNt)  char(T.Nt)   ] )' );
T.sub_nt = categorical( CharMat2CellArray([  char(T.NtSubstitute) char(T.Nt)   ] )' );
T.sub_prentntposnt = categorical( CharMat2CellArray([ char(T.NtSubstitute)   char(T.PreNt)  char(T.Nt)  char(T.PostNt) ] )' );

%% %% calculate zscores for Nt alone, prent and prent-postnt
% AM: Add mode, inferred std and z-score for each sub_prentnt
T.sub_prentnt_m_mode = NaN(height(T) , 1);
%T.sub_prentnt_s_mode = NaN(height(T) , 1);
T.sub_prentnt_z_mode = NaN(height(T) , 1);

unq_sub_prentnt = unique(T.sub_prentnt(T.IS_REF_ALLELE == 0));
for I = 1:numel(unq_sub_prentnt)
    idx_this_sub_prentnt = ismember(T.sub_prentnt , unq_sub_prentnt(I)) ;
    data_temp_sub_prentnt = T.AlleleFrequency(idx_this_sub_prentnt) ;
    m_mode = modefit(data_temp_sub_prentnt , 0 , linspace(0,prctile(T.AlleleFrequency(~T.IS_REF_ALLELE),99.9),1e3) );
    
    % use the first half (before mode) to infer distribution, symmetrical
    % around mode
    allele_freq_left_half = data_temp_sub_prentnt(data_temp_sub_prentnt <= m_mode);
	data_inferred = [allele_freq_left_half ; 2*m_mode - allele_freq_left_half];  
	s_mode = nanstd(data_inferred);
    T.sub_prentnt_m_mode(idx_this_sub_prentnt) = m_mode ; 
%    T.sub_prentnt_s_mode(idx_this_sub_prentnt) = s_mode ; 
    T.sub_prentnt_z_mode(idx_this_sub_prentnt) = (data_temp_sub_prentnt - m_mode) ./ s_mode ;
end

% for pre-nt-post tri-mer
T.sub_prentntpostnt_m_mode = NaN(height(T) , 1);
%T.sub_prentntpostnt_s_mode = NaN(height(T) , 1);
T.sub_prentntpostnt_z_mode = NaN(height(T) , 1);

unq_sub_prentntpostnt = unique(T.sub_prentntposnt(T.IS_REF_ALLELE == 0));
for I = 1:numel(unq_sub_prentntpostnt)
    idx = ismember(T.sub_prentntposnt , unq_sub_prentntpostnt(I)) ;
    data_tmp = T.AlleleFrequency(idx) ;
    m_mode = modefit(data_tmp , 0 , linspace(0,prctile(T.AlleleFrequency(~T.IS_REF_ALLELE),99.9),1e3) );

    allele_freq_left_half = data_tmp(data_tmp <= m_mode);
	data_inferred = [allele_freq_left_half ; 2*m_mode - allele_freq_left_half];  
	s_mode = nanstd(data_inferred);
    T.sub_prentntpostnt_m_mode(idx) = m_mode ; 
%    T.sub_prentntpostnt_s_mode(idx) = s_mode ; 
    T.sub_prentntpostnt_z_mode(idx) = (data_tmp - m_mode) ./ s_mode ;
end
%% add filename
T.filename = repmat( {pileupCountAllAlleles_file_name} , height(T) , 1);
T.filename = categorical(T.filename) ; 

%% shrink the table so that it takes up less space. try to keep under the 50MB github limit
T.coverage_total = [] ;
T.coverage_C = uint32( T.coverage_C ) ; 
T.coverage_W = uint32( T.coverage_W ) ; 
T.AlleleCount = uint32( T.AlleleCount ) ; 
T.Position = uint32( T.Position ) ; 

end