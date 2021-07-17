%--------------------------------------------------------------------------
% Code for the super-resolution reconstruction of video sequences
% It contains the algorithms detailed in the following publication:
% 
% 
%--------------------------------------------------------------------------


clear all
close all
close all hidden
clc


% =========================================================================

% Absolute file path of a .mat file containing the LR video in an array
% named "frames_LR" of dimension nr_lr * nc_lr * T
file_LR_vid = 'framesTeste.mat';

% upsampling factor
upsampling_fact = 2;




% =========================================================================

addpath(genpath('iterative_srr_algorithms'));
addpath(genpath('utils_funcs'));


load(file_LR_vid)

if ~exist('frames_LR')
    error('Could not load "frames_LR" from file!!!')
end

num_frames2 = size(frames_LR,3);
nr_lr = size(frames_LR,1);
nc_lr = size(frames_LR,2);




% pass "fname_rec_vid" as an additional parameter to save the video
% choose one of the SRR algorithms below:

% TSR-LMS
[frames_SR,time_reg,time_alg]=outer_super_resolve_TSR_LMS(file_LR_vid,num_frames2, nr_lr, nc_lr, upsampling_fact);
% LTSR-LMS
% [frames_SR,time_reg,time_alg]=outer_super_resolve_LTSR_LMS(file_LR_vid,num_frames2, nr_lr, nc_lr, upsampling_fact);


disp('The super-resolved sequence is in the variable named "frames_SR".')





