% created on 2017-12-19
% Goal: do STA analysis for various cells
close all;
clc;clear;
set(0,'DefaultFigureWindowStyle','docked');

exp_id = '2015_02_19';
cell_id = 'adch_62f';

data_dir = fullfile('C:\Nima\Data\',exp_id,'\');
work_dir = fullfile('.\results\',exp_id,'\');
config_file = fullfile(data_dir,'analysis_config.ini');

if ~exist(work_dir,'dir'), mkdir(work_dir); end

p = ini2struct(config_file); 

p.exp_id = exp_id;
p.cell_id = cell_id;
p.work_dir = work_dir;
p.data_dir = data_dir;

[Stimulus_all_trials,Stimulus_Matrix,V_excite,V_inhibit,Cov_Matrix,Cov_Matrix_raw_stimulus,Stim_all_trials] = STA(p);
