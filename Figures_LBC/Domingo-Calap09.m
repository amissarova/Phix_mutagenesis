%% load data
PDIR = '~/Develop/Phix_mutagenesis/' ; 
load([PDIR '/Data/T.mat']);
%% remove unneeded rows
idx = ~T.IS_REF_ALLELE & T.NtSubstitute~='a'  & T.NtSubstitute~='c'   & T.NtSubstitute~='t'   & T.NtSubstitute~='g' ;
T = T(idx,:);

%% load Rafa's data
DC = readtable([PDIR 'ExternalData/Domingo-Calap09/S2.csv'],'FileType','text');

DC.Nt = categorical(string(upper(cellfun(@(X)X(1),DC.Mutation))))  ;
DC.NtSubstitute = categorical(string(upper(cellfun(@(X)X(end),DC.Mutation))))  ;
DC.Position = cellfun(@(X)str2double(X(2:end-1)),DC.Mutation) ;

DC = innerjoin( DC, T , 'Key', {'Position' 'Nt' 'NtSubstitute'}) ; 

%%
DC.zGT2 = DC.sub_prentnt_z_mode > 2 ;
DC.zGT3 = DC.sub_prentnt_z_mode > 3 ;
DC.zGT1 = DC.sub_prentnt_z_mode > 1 ;

G = grpstats( DC ,{'Mutation' 'AminoAcidSubstitution' 'Gene'} ,{'nanmean' 'nanmedian' },'DataVars',{'RelativeFitnessEffect' 'sub_prentnt_z_mode' 'sub_prentntpostnt_z_mode' 'AlleleFrequency' 'zGT3' 'zGT2' 'zGT1'}) ;
G = sortrows( G ,{ 'nanmean_RelativeFitnessEffect' 'nanmean_sub_prentnt_z_mode'});
G = G(G.nanmean_RelativeFitnessEffect > -0.5 ,:);
G = G(~isnan(G.nanmean_sub_prentnt_z_mode) & ~isinf(G.nanmean_sub_prentnt_z_mode)   , :);
figure;
X = G.nanmean_RelativeFitnessEffect  ; 
Y = G.nanmean_sub_prentnt_z_mode ; 
[c,p] = corr(X,Y);
gscatter( X , Y , G.Gene , [] , 'vop+');
xlabel('Fitness (Domingo-Calap ''09)')
ylabel('Z score')
title(sprintf('c = %0.02f p=%0.04f' , c , p ));
