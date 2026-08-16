function main_cameraparams(raw_fmri_image_path, raw_fmri_image_path_bgremoved, output_dir, ...
                           motion_params_file, displacement_file, slice_num_ordering_file, ...
                           recon_noscrub_filename_prefix, recon_scrub_filename_prefix, ...
                           thresh_forscrub, sms_fac)
%% Example Usage Via Command Line:
% matlab -batch \
% "main_cameraparams( \
% '/lab-share/Neuro-Cohen-e2/Groups/IRB-P00049401/intravolume_motion_pipeline/participant_data/func-bold_task-NFB2_20250827181314_24.nii.gz', \
% '/lab-share/Neuro-Cohen-e2/Groups/IRB-P00049401/intravolume_motion_pipeline/participant_data/func-bold_task-NFB2_20250827181314_24_bgremoved.nii.gz', \
% '/lab-share/Neuro-Cohen-e2/Groups/IRB-P00049401/intravolume_motion_pipeline/test', \
% '/lab-share/Neuro-Cohen-e2/Groups/IRB-P00049401/intravolume_motion_pipeline/participant_data/radian-parameters.txt', \
% '/lab-share/Neuro-Cohen-e2/Groups/IRB-P00049401/intravolume_motion_pipeline/participant_data/displacements.txt', \
% '/lab-share/Neuro-Cohen-e2/Groups/IRB-P00049401/intravolume_motion_pipeline/participant_data/slice_num_ordering_sms-fac-4.txt', \
% 'no-scrubbing_motion-corrected_func_image', \
% 'scrubbed_and_motion-corrected_func_image', \
% 0.4157, \
% 4 \
% )"
%%

close all;
clc;
set(0,'DefaultFigureWindowStyle','docked');
%%  Add all paths
addpath('../direct-liftandunlift-codes');
addpath('../operators/');
addpath('../');
%% create folders to save background volumes and motion params

opt.output_dir = output_dir;
[~, ~, ~] = mkdir(output_dir);

if ischar(thresh_forscrub) || isstring(thresh_forscrub)
    thresh_forscrub = str2double(thresh_forscrub);
end

if ischar(sms_fac) || isstring(sms_fac)
    sms_fac = str2double(sms_fac);
end

Y = py.SimpleITK.ReadImage(raw_fmri_image_path,py.SimpleITK.sitkFloat64);
origin_4d = Y.GetOrigin();
spacing_4d = Y.GetSpacing();
direction_4d=  Y.GetDirection();
%Iorig = numpytomatlab(py.SimpleITK.GetArrayFromImage(Y));
Iorig =  permute(double(py.SimpleITK.GetArrayFromImage(Y)),[3,4,2,1]);
[n1,n2,nsl,nvs] = size(Iorig);
opt.vol_start = 1; % Adjust this parameter to discard initial frames, vol_start =  No of FramesTodiscard + 1
nv = nvs-opt.vol_start+1;
if opt.vol_start ~=1
    Iorig_tr = Iorig(:,:,:,opt.vol_start:nvs);
    Iorig_tr_np = py.numpy.asarray(permute(Iorig_tr,[4,3,1,2]));
    Y = py.pyfuncs_forsltovol.numpy4Dtositk(Iorig_tr_np,origin_4d,direction_4d,spacing_4d,int32(nv));
end
opt.I1 = Iorig(:,:,:,opt.vol_start);

%Iorig_tr = Iorig(:,:,:,vol_start:nvs);
%Inp =  matlabtonumpy(I);
%Inp = py.numpy.asarray(permute(Iorig_tr,[4,3,1,2]));
%Y =  py.vvr_regtofirstvolofmo.numpy4Dtositk(Inp,origin_4d,direction_4d,spacing_4d,int32(nv));
opt.nv = nv;
nslbysmsfac = nsl/sms_fac;
opt.nslbysmsfac = nslbysmsfac;
opt.sms_fac = sms_fac;

opt.slice_acq_order = load(slice_num_ordering_file)';
opt.slice_acq_order = opt.slice_acq_order + 1;
assert(numel(opt.slice_acq_order) == nsl, ...
    'Slice acquisition order must contain exactly %d entries.', nsl);

assert(isequal(sort(opt.slice_acq_order), 1:nsl), ...
    'Slice acquisition order must contain every slice exactly once.');

opt.params_ind = py.numpy.int32(py.range(int32(nslbysmsfac)));
const = sms_fac*ones(1,nslbysmsfac);
sl_acq_order_cell_py = mat2cell(opt.slice_acq_order-1,[1],const); % python convention.
%ilacq_cell_mat = mat2cell(opt.slice_acq_order,[1],const);  % matlab convention.
slice_info = py.dict(pyargs());
for i = 1:nslbysmsfac
    slice_info{i-1} =  matlabtonumpy(int32(sl_acq_order_cell_py{i}));% according to python convention.
end
opt.slice_info = slice_info;

%% Compute framewise displacement
params_temp = load(motion_params_file);
[r,~] = size(params_temp);
params_temp = params_temp(opt.vol_start:r,:);

% Load parameters from radian_parameters.txt: to be used for slice-level motion correction in reconstruct_timeseries_noscrub_forsltovol
opt.params = py.numpy.asarray(params_temp);

% Load displacements from displacements.txt: to be used to determine which volumes to be scrubbed in MainFunctionForInterpScrubbedData
SWD = load(displacement_file);

%% Create SWD plot without displaying a GUI
set(groot, 'defaultFigureWindowStyle', 'normal');
set(groot, 'defaultFigureVisible', 'off');

fig = figure('Visible', 'off', 'WindowStyle', 'normal');

plot(SWD);
hold on;

y = thresh_forscrub * ones(length(SWD), 1);
plot(y);

xlabel('Aquisition');
ylabel('Intravolume-Motion');
title('Intravolume-Motion');

exportgraphics(fig, fullfile(output_dir, 'SWD.png'), 'Resolution', 300);

close(fig);

disp("SWD Plot Saved to:")
disp(fullfile(output_dir, 'SWD.png'))
%%
opt.n1 = n1;
opt.n2 = n2;
opt.nsl = nsl;
%% fmri time series without background
%raw_fmri_image_path_bgremoved = strcat(data_path,'_bgremoved.nii.gz');
Y_nobg_img = py.SimpleITK.ReadImage(raw_fmri_image_path_bgremoved,py.SimpleITK.sitkFloat64);
I_nobg = permute(double(py.SimpleITK.GetArrayFromImage(Y_nobg_img)),[3,4,2,1]);
%I_nobg = numpytomatlab(py.SimpleITK.GetArrayFromImage(Ynobg_img));
I_nobg = I_nobg(:,:,:,opt.vol_start:nvs);
X_nobg = reshape(I_nobg,[n1*n2*nsl,nv]);
opt.ind_bg = ~any(X_nobg,2); % index value of 1 corresponds to background pixel.
clear I_nobg X_nobg Y_nobg_img
%% set parameters
opt.mu = 25; % regularization parameter -  From my experiments, I have observed the algorithm to be robust with changing mu.
%% For now, don't change the value of the pair (beta, beta_fac). Generate the results and we can decide if we want to try a different pair.
opt.beta = 7.5;%0.25;%7.5;% initial value of beta
opt.beta_fac = 1.1;%1.2; %1.1;% variable to increment beta every iteration
%%
opt.ft = 50;%75;% filter size. Determines the size of Hankel matrix formed at every voxel.
opt.overall_maxIter = 20;%15; % Maximum number of iterations for the algorithm
opt.maxIter = 6;% Maximum number of iterations for the x-subproblem
opt.Njobs = -1;% Uses all cores when the python codes corresponding to the z-subproblem are executed
opt.timepoint1 = 116819;% For display purposes. time series corresponding to it is displayed
opt.timepoint2 = 104879;% For display purposes. time series corresponding to it is displayed
opt.ts_start  = 1; % Starting point for downsampling the reconstructed time series
opt.tolerance = 1e-4; % Stopping point for X subproblem
opt.eta = 1.1;% parameter for IRLS algorithm. Usually, eta > 1 and < 1.5opt.I1 = Iorig(:,:,:,opt.vol_start);  % for initialization purposes.
opt.p = 0.1;% Schatten p-norm value 0 <= p <=1


[~,X_img,opt] = reconstruct_timeseries_noscrub_forsltovol(Y,opt);


%fname = strcat(data_path,'recon_',g,'.nii.gz');
%fname = strcat(data_path,'recon_',g,'_beta0.2_fac1.2_iter20_segnooverlap.nii.gz');
fpath = fullfile(output_dir, [recon_noscrub_filename_prefix '.nii.gz']);
py.SimpleITK.WriteImage(X_img, fpath);  
disp("Reconstructed fMRI Image (Pre-Scrubbing) at: ")
disp(fpath)

%fname = strcat(data_path,'recon_noscrub_betai nit',num2str(opt.beta),'_',g,'.nii.gz');
    %py.SimpleITK.WriteImage(X_img, fname);
%X_img = py.SimpleITK.ReadImage(fname,py.SimpleITK.sitkFloat64);
%% The downsampled output is the input to the second stage, where volumes with FD > threshold are marked for scrubbing.
%% The missing data is then interpolated using the matrix completion algorithm.
%% Interpolate scrubbed data.


%% Scrub volumes if ANY acquisition-to-acquisition movement exceeds threshold

nGroups = nslbysmsfac;   %ex. 15
nVolumes = nv;           %ex. 283

assert(length(SWD) == nGroups*nVolumes - 1, ...
    'Unexpected number of displacement values.');

max_displacement_per_volume = zeros(nVolumes,1);

for v = 1:nVolumes
    start_idx = (v-1)*nGroups + 1;
    end_idx = min(v*nGroups, length(SWD));

    max_displacement_per_volume(v) = max(SWD(start_idx:end_idx));
end

opt.lower_thresh = -1;
opt.higher_thresh = 2;
opt.tolerance = 1e-6;
opt.maxIter = 45;
opt.thresh_forscrub = thresh_forscrub;
opt.SWD = max_displacement_per_volume;
opt.beta = 1; %% Don't change.

fig = figure('Visible','off');
plot(max_displacement_per_volume);
hold on;

y = thresh_forscrub * ones(1,length(max_displacement_per_volume));
plot(y);

xlabel('Volume');
ylabel('Maximum intra-volume displacement');
title('Volume-wise Motion / Scrubbing Threshold');

filename = fullfile(output_dir, 'max-swd-per-volume.png');
exportgraphics(fig, filename, 'Resolution', 300);
close(fig);

disp("Scrubbing Plot Saved to:")
disp(fullfile(filename))

[vol_i,~,~] = find(max_displacement_per_volume>=thresh_forscrub);

if ~isempty(vol_i)
    opt.filt_siz = round(nv/2);   %% Have set the filter size to length of series/2. If the results are not good, run the algorithm for a different value.
    for count = 1:length(opt.filt_siz)
        opt.ft = opt.filt_siz(count);
        [Xs,vol_rem] = MainFunctionForInterpScrubbedData(X_img,opt);
        Is = reshape(Xs,n1,n2,nsl,nv);
        X_np = matlabtonumpy(Is);
        X_img = py.pyfuncs_sv_parallelized_ss_volumelevel.numpy4Dtositk(X_np,origin_4d,direction_4d,spacing_4d,int32(nv));
        filename = fullfile(output_dir, [recon_scrub_filename_prefix '.nii.gz']);
        py.SimpleITK.WriteImage(X_img, filename);
        disp("Reconstructed fMRI Image (Post-Scrubbing) at: ")
        disp(strcat(filename))
    end
    py.numpy.save(strcat(output_dir,'vol_rem_afterscrubbing'),matlabtonumpy(vol_rem));
end
disp(strcat(g, '    done'));
