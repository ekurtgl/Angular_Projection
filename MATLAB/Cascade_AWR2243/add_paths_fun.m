function [] = add_paths_fun()
%     homeDir = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO'); % for windows
    homeDir = '/mnt/HDD04/Projection_data/Scripts/Cascade_AWR2243/MatlabExamples/4chip_cascade_MIMO_example';
    addpath(genpath([homeDir,'/modules']));
    addpath(genpath([homeDir,'/main']));
    addpath([homeDir,'/utils/math']);
    addpath([homeDir,'/utils/dataParse']);
    addpath([homeDir,'/utils/disp']);
    addpath([homeDir,'/utils/cascade_json_parser']);
end