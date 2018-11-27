%% load data 
addpath(genpath('~/Develop/Phix_mutagenesis/'));
%load Data/T.mat
%%
T = GenerateSaveTableFromAll_pileupCountAlleles_tabfiles(); 
%%
load ../Data/Stat
D.Nt = categorical(D.Nt);
D.NtSubstitute = categorical(D.NtSubstitute);
D.IN_ORF = categorical(D.IN_ORF);
D.SynonymousTotalBool = categorical(D.SynonymousTotalBool);
D.NonsenseTotalBool = categorical(D.NonsenseTotalBool);
D.MutationType = categorical(D.MutationType);
D.ConservationScore = categorical(D.ConservationScore);

vn = {'SynonymousTotalBool' 'NonsenseTotalBool' 'MutationType' 'ConservationScore' 'IN_ORF' 'Position' 'Nt' 'NtSubstitute'};
R = innerjoin( T , dataset2table(D(:,vn)) ,  'Key', {'Position' 'Nt' 'NtSubstitute'});

%% make per-file matrix of z scores
% % fn = regexprep( string(T.filename) , '^.*\/','');
% % fn = regexprep( fn , '\.tab$','');
% % T.fn = regexprep( fn , '_',' ');

T = sortrows(T,{'filename' 'Position' 'Nt' 'NtSubstitute'},'ascend');
ufn = unique(T.filename);
%% build matrix
data = T.AlleleFrequency ; 

dmat = reshape( data , [] , numel(ufn));


figure ; 
imagesc(corr(dmat,'rows','complete','type','Spearman') , [0.4 1])
colorbar
title('Allele Frequency')

%%
R.zGT3 = R.sub_prentnt_z_mode > 3 ; 

modelspec = 'zGT3 ~ SynonymousTotalBool+NonsenseTotalBool+MutationType+ConservationScore+IN_ORF';

mdl = fitglm(R,modelspec,'Distribution','binomial')

%%
figure; hold on ; 
ecdf( R.sub_prentnt_z_mode( R.IN_ORF=='1')  );
ecdf( R.sub_prentnt_z_mode( R.IN_ORF=='0') );
xlim([-2 5])
%%
R.sub_prentnt_z_mode(R.sub_prentnt_z_mode>5)=5;
R.data = R.sub_prentnt_z_mode ; 
R.data = R.data>3 ;

figure;  
%G = grpstats( R , { 'fn' 'SynonymousTotalBool'},{'mean' 'median' 'std' 'sem'},'DataVars','data');
G = grpstats( R , { 'fn' 'ConservationScore'},{'mean' 'median' 'std' 'sem'},'DataVars','data');

gscatter( G.mean_data(5:6:end) , G.mean_data(2:6:end) , G.fn(2:6:end))
xlabel( char(G.ConservationScore(5)) )
ylabel( char(G.ConservationScore(2)) )
%xlim([-0.5 0.75])
%ylim(xlim)
line(xlim,xlim)