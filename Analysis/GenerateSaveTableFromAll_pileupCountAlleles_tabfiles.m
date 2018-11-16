function T = GenerateSaveTableFromAll_pileupCountAlleles_tabfiles(  ) 
%% load all Illumina pileupCountAllAlleles tab files, build and save the dataset
% 
%
% 

PROJECT_DIR = '~/Develop/Phix_mutagenesis/';
addpath_recurse( PROJECT_DIR );
FILEREGEXP = [ PROJECT_DIR '/Illumina_Data/*/*.tab' ];
file_list = dir(FILEREGEXP);
file_list_cell_array = strcat({file_list.folder} ,filesep,{file_list.name}) ; 

T = table();
for fileI = 1:numel(file_list_cell_array)
    fprintf('%d/%d\t%s\n' , fileI , numel(file_list_cell_array) ,  file_list_cell_array{fileI}) ;
    T = vertcat( T , reshape_pileupCountAllAlleles(file_list_cell_array{fileI}) ) ;
end

save( [PROJECT_DIR 'Data/T.mat'] , 'T' );

end


