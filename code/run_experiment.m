% Copyright (c) 2026 Caiming Tan
% SPDX-License-Identifier: MIT

%RUN_EXPERIMENT Reproduce the three-position mixed-accuracy IMU comparison.
clc
clearvars

project_root = fileparts(fileparts(mfilename('fullpath')));
code_dir = fileparts(mfilename('fullpath'));
addpath(code_dir);

data_dir = fullfile(project_root, 'data', 'experimental', '2025-09-04', 'raw');
zimu_file = fullfile(data_dir, 'zimu070_20250904_203714.txt');
mems_file = fullfile(data_dir, 'er3mg01_20250904_203718.txt');
zimu_data = load(zimu_file);
mems_data = load(mems_file);

earth_rate = 7.2921158e-5;
latitude = 32.113350 * pi / 180;
gravity = 9.7803253 * (1 + 0.0053022 * sin(latitude)^2 ...
    - 0.0000058 * sin(2 * latitude)^2);
vector_angle = pi / 2 - latitude;
reference_p = [0; 0; 1];
reference_q = [0; cos(latitude); sin(latitude)];
body_axis_transform = [0,1,0;-1,0,0;0,0,1];
time_origin_ms = min(zimu_data(1,1), mems_data(1,1));
zimu_elapsed_ms = zimu_data(:,1) - time_origin_ms;
mems_elapsed_ms = mems_data(:,1) - time_origin_ms;
position_windows_ms = [2e4,6e4; 12e4,16e4; 24e4,28e4];
position_names = ["Position 1"; "Position 2"; "Position 3"];

unknown_index_i = 3;
observed_index_j = 1;
sigma_a_i1 = 1e-3;
sigma_a_i2 = 1e-3;
sigma_b_j = 2e-3;
selection_threshold = 3 * sigma_b_j;
results = zeros(size(position_windows_ms,1), 8);

for position_index = 1:size(position_windows_ms,1)
    window = position_windows_ms(position_index,:);
    zimu_mask = zimu_elapsed_ms >= window(1) & zimu_elapsed_ms <= window(2);
    mems_mask = mems_elapsed_ms >= window(1) & mems_elapsed_ms <= window(2);
    if ~any(zimu_mask) || ~any(mems_mask)
        error('No samples found in the time window for %s.', ...
            position_names(position_index));
    end

    zimu_mean = mean(zimu_data(zimu_mask,2:7),1);
    zimu_gyro = zimu_mean(1:3)' * pi / 180;
    accelerometer = zimu_mean(4:6)';
    mems_mean = mean(mems_data(mems_mask,2:4),1);
    mems_gyro = [-mems_mean(2), -mems_mean(1), mems_mean(3)]' * pi / 180;

    observed_a = body_axis_transform * (accelerometer / gravity);
    reference_b_fog = body_axis_transform * (zimu_gyro / earth_rate);
    observed_b_mems = body_axis_transform * (mems_gyro / earth_rate);

    ground_truth_dcm_analysis = estimateTriadAttitude( ...
        observed_a,reference_b_fog,reference_p,reference_q);
    ground_truth_dcm = body_axis_transform' * ground_truth_dcm_analysis;
    ground_truth_angles = cnbToEulerZxy(ground_truth_dcm') * 180 / pi;

    a_i1 = observed_a(cyclicIndex(unknown_index_i,1));
    a_i2 = observed_a(cyclicIndex(unknown_index_i,2));
    b_j = reference_b_fog(observed_index_j);
    [a_candidates,b_candidates] = constructMiadCandidates( ...
        unknown_index_i,observed_index_j,a_i1,a_i2,b_j,vector_angle, ...
        sigma_a_i1,sigma_a_i2,sigma_b_j);

    a_i = observed_a(unknown_index_i);
    b_j1 = observed_b_mems(cyclicIndex(observed_index_j,1));
    b_j2 = observed_b_mems(cyclicIndex(observed_index_j,2));
    [miad_a,miad_b] = selectCandidateWithFullInformation( ...
        unknown_index_i,observed_index_j,a_candidates,b_candidates,a_i,b_j1,b_j2);
    mixed_b = observed_b_mems;
    mixed_b(observed_index_j) = b_j;

    miad_dcm_analysis = estimateTriadAttitude( ...
        miad_a,miad_b,reference_p,reference_q);
    triad_dcm_analysis = estimateTriadAttitude( ...
        observed_a,mixed_b,reference_p,reference_q);
    weights = [sigma_b_j^2, sigma_a_i1^2];
    weights = weights / sum(weights);
    quaternion = solveQuestAttitude(observed_a,mixed_b,reference_p,reference_q, ...
        weights(1),weights(2));
    quest_cnb_analysis = quaternionToCnb(quaternion)';
    quest_dcm_analysis = quest_cnb_analysis';

    miad_dcm = body_axis_transform' * miad_dcm_analysis;
    triad_dcm = body_axis_transform' * triad_dcm_analysis;
    quest_dcm = body_axis_transform' * quest_dcm_analysis;
    miad_angles = cnbToEulerZxy(miad_dcm') * 180 / pi;
    triad_angles = cnbToEulerZxy(triad_dcm') * 180 / pi;
    quest_angles = cnbToEulerZxy(quest_dcm') * 180 / pi;

    gt_az = round(ground_truth_angles(3)*100)/100;
    tr_az = round(triad_angles(3)*100)/100;
    qs_az = round(quest_angles(3)*100)/100;
    md_az = round(miad_angles(3)*100)/100;
    azimuth_error = [tr_az - gt_az,  qs_az - gt_az,  md_az - gt_az];

    results(position_index,:) = [position_index, gt_az, tr_az, qs_az, md_az, azimuth_error];
end

result_table = array2table(results, 'VariableNames', ...
    {'Position','GroundTruthDeg','TriadDeg','QuestDeg','MiadDeg', ...
    'TriadErrorDeg','QuestErrorDeg','MiadErrorDeg'});
result_table.Position = position_names;
disp(result_table);

results_dir = fullfile(project_root, 'results', 'experiment');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
writetable(result_table,fullfile(results_dir,'experiment_results.csv'));
save(fullfile(results_dir,'experiment_results.mat'), ...
    'result_table','position_windows_ms','unknown_index_i','observed_index_j');
