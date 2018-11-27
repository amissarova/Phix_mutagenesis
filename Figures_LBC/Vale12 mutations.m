%% load data
PDIR = '~/Develop/Phix_mutagenesis/' ; 
load([PDIR '/Data/T.mat']);
%% remove unneeded rows
idx = ~T.IS_REF_ALLELE & T.NtSubstitute~='a'  & T.NtSubstitute~='c'   & T.NtSubstitute~='t'   & T.NtSubstitute~='g' ;
T = T(idx,:);

%% load Rafa's data
VALE1 = readtable([PDIR 'ExternalData/Vale12/Table1.tsv'],'FileType','text');
VALE2 = readtable([PDIR 'ExternalData/Vale12/Figure1.tab'],'FileType','text');
VALE = join(VALE1,VALE2,'Key','MutantID');
clear 'VALE1' 'VALE2'

VALE.Nt = categorical(string(upper(cellfun(@(X)X(1),VALE.Mutation))))  ;
VALE.NtSubstitute = categorical(string(upper(cellfun(@(X)X(end),VALE.Mutation))))  ;
VALE.Position = cellfun(@(X)str2double(X(2:end-1)),VALE.Mutation) ;

VALE = innerjoin( VALE, T , 'Key', {'Position' 'Nt' 'NtSubstitute'}) ; 
VALE.filename = char(VALE.filename);


% any relation between our data and theirs? 
VALE.zGT2 = VALE.sub_prentnt_z_mode > 2 ;
VALE.zGT3 = VALE.sub_prentnt_z_mode > 3 ;
VALE.zGT1 = VALE.sub_prentnt_z_mode > 1 ;

%VALE.sub_prentnt_z_mode(VALE.sub_prentnt_z_mode > 10)=10 ; 
%VALE.sub_prentntpostnt_z_mode(VALE.sub_prentntpostnt_z_mode > 10)=10 ; 
%VALE.sub_prentnt_z_mode(VALE.sub_prentnt_z_mode < -10)=-10 ; 
%VALE.sub_prentntpostnt_z_mode(VALE.sub_prentntpostnt_z_mode < -10)=-10 ;

idx_to_drop = regexpcmp(VALE.filename,'extendedFrags') | regexpcmp(VALE.filename,'NovaSeq') ; 
%idx_to_drop = VALE.coverage_W < 6e5 ; 
G = grpstats(VALE(~idx_to_drop,:),{'Mutation' 'MutantID' 'AminoAcidSubstitution'} ,{'nanmean' 'nanmedian' },'DataVars',{'Fitness' 'sub_prentnt_z_mode' 'sub_prentntpostnt_z_mode' 'AlleleFrequency' 'zGT3' 'zGT2' 'zGT1'}) ;
G = sortrows(G,'MutantID');

figure;
X = G.nanmedian_Fitness  ; 
Y = G.nanmean_zGT2 ;
[c,p] = corr(X,Y);
gscatter( X , Y , G.MutantID , [] , 'ov+p.' )
xlabel('Fitness (Vale ''12)')
ylabel('Z score')
title(sprintf('c = %0.02f p=%0.04f' , c , p ));
