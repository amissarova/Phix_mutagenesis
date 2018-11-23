%% load data
PDIR = '~/Develop/Phix_mutagenesis/' ; 
load([PDIR '/Data/T.mat']);

% set labels and groups for boxplotting allele freqs
idx = ~T.IS_REF_ALLELE & T.NtSubstitute~='a'  & T.NtSubstitute~='c'   & T.NtSubstitute~='t'   & T.NtSubstitute~='g' ;
T = T(idx,:);
sym = strcat(string(T.Nt) , '->' , string(T.NtSubstitute) );
sym2 = strcat( string(T.PreNt) ,  string(T.Nt) , '>' , string(T.NtSubstitute) );

%% raw allele freqs
figure;
boxplot(log10(T.AlleleFrequency ),{T.Nt T.NtSubstitute T.PreNt},'Symbol','','Labels',(sym2),'notch','on');
ylim([min(ylim) -3])
ylabel('log10( allele frequency )')
xlabel('Substitution')
%% one example
% AC->AG & CC->CG
idxA = T.Nt == 'C'  & T.NtSubstitute == 'G' & T.PreNt == 'A'  ;
idxC = T.Nt == 'C'  & T.NtSubstitute == 'G' & T.PreNt == 'C'  ;

fh = figure('units','centimeters','position',[5 5 9 8]);
hold on ; 
data =  T.AlleleFrequency ; 
data(data>1e-4) = 1e-4 ; 
histogram( ( data( idxA ) ) , 50 , 'Normalization','Probability')
histogram( ( data( idxC ) ) , 50 , 'Normalization','Probability')
ylabel('Fraction of positions & experiments')
xlabel('Allele frequency')
legend( {'aC -> aG' 'cC -> cG'})

fh = figure('units','centimeters','position',[5 5 9 8]);
hold on ; 
data =  T.sub_prentnt_z_mode ;
xl = linspace(-5,5,50);
data(data>5) = 5 ; 
histogram( ( data( idxA ) ) , xl , 'Normalization','Probability')
histogram( ( data( idxC ) ) , xl , 'Normalization','Probability')
ylabel('Fraction of positions & experiments')
xlabel('"z-score" a.f. using left half from mode')
legend( {'aC -> aG' 'cC -> cG'})
xlim([-3 5])
%% raw z-score
figure;
boxplot( T.sub_prentnt_z_mode ,{T.Nt T.NtSubstitute T.PreNt},'Symbol','','Labels',(sym2),'notch','on');
ylim([-2 10 ])
ylabel('"z-score" a.f. using left half from mode')
xlabel('Substitution')

%% across experiments;
T = sortrows( T, {'filename' , 'Position' ,'NtSubstitute'});
data = reshape( T.AlleleFrequency , max(T.Position)*3 , []);
c = corr(data,'type','Spearman','rows','complete');

fh = figure('units','centimeters','position',[5 5 17 15]);
imagesc(c);
colorbar;
xlabel('experiment')
ylabel('experiment')
title('Spearman correlation of z-scores across experiments')

%%
fn = string(T.filename) ;
fn = regexprep( fn(1:(max(T.Position)*3):end)  , '.*\/Illumina_Data\/' , '')
fn = regexprep( fn , '\..*','')
fn = regexprep( fn , 'ileupCoun.*','')

%%
idx = ~cellfun(@isempty,(regexp(fn,'_[lL]'))) ;
c = corr(data,'type','Spearman','rows','complete');
c = c(idx,idx);
fh = figure('units','centimeters','position',[5 5 17 15]);
imagesc(c,[0.5 1]);
colorbar;
xlabel('experiment')
ylabel('experiment')
title('Spearman correlation of Allele Freqs across experiments')
