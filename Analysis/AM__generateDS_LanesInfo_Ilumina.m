%% 1. Pull the template DS_mutations_info w/ info regarding every possible substitution
cd ~/Develop/Phix_mutagenesis/Illumina_Data
addpath(genpath('~/Develop/matlab'));
load('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat'); 
D = DS;
Illumina_data = dir('~/Develop/Phix_mutagenesis/Illumina_Data'); 
for I = 1:length(Illumina_data)
    if isfolder(Illumina_data(I).name) & ~strcmp(Illumina_data(I).name , '.') & ~strcmp(Illumina_data(I).name , '..')
        D = addPhixMutationIlluminaDataset(D , Illumina_data(I).name);
    end
end
%save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 