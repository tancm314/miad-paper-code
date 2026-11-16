% Copyright (c) 2026 Caiming Tan
% SPDX-License-Identifier: MIT

%RUN_SIMULATION Reproduce the Monte Carlo results and RMSE figure.
clc
close all
clearvars

project_root = fileparts(fileparts(mfilename('fullpath')));

rng(2, 'twister');

sigma_a_i1 = 1e-3;
sigma_a_i2 = 1e-3;
sigma_b_j = 1e-3;
sigma_a_i = 1e-2;
sigma_b_j1 = 1e-2;
sigma_b_j2 = 1e-2;
threshold = 3 * sigma_b_j2;
vector_angle = 45 * pi / 180;
num_test_cases = 50;
num_monte_carlo_samples = 1000;
results_dir = fullfile(project_root, 'results', 'simulation');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

M = num_test_cases;
save_RMSE_element_meas=zeros(M,3);
save_RMSE_element_theory=zeros(M,3);
save_RMSE_phi_meas=zeros(M,3);
save_RMSE_phi_triad_meas=zeros(M,3);
save_RMSE_phi_quest_meas=zeros(M,3);
save_RMSE_phi_theory=zeros(M,3);
save_RMSE_phi_triad_theory=zeros(M,3);
save_RMSE_phi_quest_theory=zeros(M,3);
save_RMSE_dis_meas=zeros(M,1);
save_RMSE_dis_triad_meas=zeros(M,1);
save_RMSE_dis_quest_meas=zeros(M,1);
save_RMSE_dis_theory=zeros(M,1);
save_theory=zeros(M,3);
save_atti_true=zeros(M,3);
save_theta_true=zeros(M,1);
save_N_count_crit1=zeros(M,4);
save_N_count_crit2=zeros(M,4);
save_mean_dis1=zeros(M,2);
save_mean_dis2=zeros(M,2);
save_i_j_a_b=zeros(M,9);
save_exc=zeros(M,21);
save_norm_violation_count=zeros(M,1);
save_delta_revision_count=zeros(M,2);
for k = 1:M
    N = num_monte_carlo_samples;
    i = randperm(3,1);
    j = randperm(3,1);
    p = [1;0;0];
    q = [cos(vector_angle);sin(vector_angle);0];
    [observed_a_samples,observed_b_samples,a_true,b_true,true_attitude] = ...
        generateVectorMeasurements(i,j,sigma_a_i1,sigma_a_i2,sigma_b_j, ...
        sigma_a_i,sigma_b_j1,sigma_b_j2,p,q,N);
    Cnb_true=eulerZxyToCnb(true_attitude);
    Cbn_true=Cnb_true';

    s_element=zeros(N,3);
    s_phi=zeros(N,3);
    s_phi_triad=zeros(N,3);
    s_phi_quest=zeros(N,3);
    s_dis=zeros(N,1);
    s_dis_triad=zeros(N,1);
    s_dis_quest=zeros(N,1);
    s_crit1=zeros(N,3);
    s_crit2=zeros(N,3);
    s_exc=zeros(N,3);
    s_exc_value=zeros(N,3);
    s_epsilon=zeros(N,1);
    s_epsilon2=zeros(N,1);
    for n = 1:N
        a = observed_a_samples(n,:)';
        b = observed_b_samples(n,:)';
        a_i1 = a(cyclicIndex(i,1));
        a_i2 = a(cyclicIndex(i,2));
        b_j = b(j);
        [a_candidates,b_candidates,exc,exc_value,epsilon,epsilon2] = ...
            constructMiadCandidates(i,j,a_i1,a_i2,b_j,vector_angle, ...
            sigma_a_i1,sigma_a_i2,sigma_b_j);

        a_i = a(i);
        sign_a_i = sign(a_i);
        b_j1 = b(cyclicIndex(j,1));
        b_j2 = b(cyclicIndex(j,2));
        [a_criterion1,b_criterion1,alarm1] = selectCandidateWithMinimalInformation( ...
            i,j,a_candidates,b_candidates,sign_a_i,b_j1,b_j2,true,threshold);
        [a_criterion2,b_criterion2,alarm2] = selectCandidateWithMinimalInformation( ...
            i,j,a_candidates,b_candidates,sign_a_i,b_j1,b_j2,false,threshold);
        [selected_a,selected_b] = selectCandidateWithFullInformation( ...
            i,j,a_candidates,b_candidates,a_i,b_j1,b_j2);

        Cbn=estimateTriadAttitude(selected_a,selected_b,p,q);
        Cbpb=Cbn*Cnb_true;
        phi=[-Cbpb(2,3);Cbpb(1,3);-Cbpb(1,2)]; 
        dis=computeAttitudeDistance(Cbn,Cbn_true,2);

        Cbn_triad=estimateTriadAttitude(a,b,p,q);
        Cbpb_triad=Cbn_triad*Cnb_true;
        phi_triad=[-Cbpb_triad(2,3);Cbpb_triad(1,3);-Cbpb_triad(1,2)];
        dis_triad=computeAttitudeDistance(Cbn_triad,Cbn_true,2);

        weights=[sigma_b_j1^2,sigma_a_i1^2];
        weights=weights/sum(weights);
        quaternion = solveQuestAttitude(a,b,p,q,weights(1),weights(2));
        Cnb_quest=quaternionToCnb(quaternion)';
        Cbn_quest=Cnb_quest';
        Cbpb_quest=Cbn_quest*Cnb_true;
        phi_quest=[-Cbpb_quest(2,3);Cbpb_quest(1,3);-Cbpb_quest(1,2)]; 
        dis_quest=computeAttitudeDistance(Cbn_quest,Cbn_true,2);

        Cbn1=estimateTriadAttitude(a_criterion1,b_criterion1,p,q);
        dis_atti1=computeAttitudeDistance(Cbn1,Cbn,2);
        Cbn2=estimateTriadAttitude(a_criterion2,b_criterion2,p,q);
        dis_atti2=computeAttitudeDistance(Cbn2,Cbn,2);
        dis_vec1=norm(extractDiscriminatorComponents( ...
            i,j,selected_a,selected_b)- ...
            extractDiscriminatorComponents(i,j,a_criterion1,b_criterion1));
        dis_vec2=norm(extractDiscriminatorComponents( ...
            i,j,selected_a,selected_b)- ...
            extractDiscriminatorComponents(i,j,a_criterion2,b_criterion2));

        s_element(n,:)=[selected_a(i)-a_true(i), ...
            selected_b(cyclicIndex(j,1))-b_true(cyclicIndex(j,1)), ...
            selected_b(cyclicIndex(j,2))-b_true(cyclicIndex(j,2))];
        s_phi(n,:)=phi;
        s_phi_triad(n,:)=phi_triad;
        s_phi_quest(n,:)=phi_quest;
        s_dis(n)=dis;
        s_dis_triad(n)=dis_triad;
        s_dis_quest(n)=dis_quest;    
        s_crit1(n,:)=[dis_vec1,dis_atti1,alarm1];
        s_crit2(n,:)=[dis_vec2,dis_atti2,alarm2];
        s_exc(n,:)=exc;
        s_exc_value(n,:)=exc_value;
        s_epsilon(n)=epsilon;
        s_epsilon2(n)=epsilon2;
    end
    P_element_meas=s_element'*s_element/N;
    P_phi_meas=s_phi'*s_phi/N;
    P_phi_triad_meas=s_phi_triad'*s_phi_triad/N;
    P_phi_quest_meas=s_phi_quest'*s_phi_quest/N;
    RMSE_dis_meas=sqrt(s_dis'*s_dis/N);
    RMSE_dis_triad_meas=sqrt(s_dis_triad'*s_dis_triad/N);
    RMSE_dis_quest_meas=sqrt(s_dis_quest'*s_dis_quest/N);

    [RMSE_element,RMSE_phi,RMSE_dis,denom,denom2]=propagateMiadUncertainty(i,j,...
        a_true,b_true,sigma_a_i1,sigma_a_i2,sigma_b_j);
    RMSE_phi_triad=estimateTriadRmse(a_true,b_true,sigma_a_i1,sigma_b_j);
    RMSE_phi_quest=estimateQuestRmse(a_true,b_true,sigma_a_i1,sigma_b_j);

    [N_count_crit1,mean_dis_atti1,mean_dis_vec1]=countMiadStates(s_crit1,1);
    [N_count_crit2,mean_dis_atti2,mean_dis_vec2]=countMiadStates(s_crit2,2);

    nonzero_epsilon = s_epsilon(s_epsilon ~= 0);
    if isempty(nonzero_epsilon)
        epsilon_min=0;
        epsilon_max=0;
    else
        epsilon_min=min(nonzero_epsilon);
        epsilon_max=max(nonzero_epsilon);
    end

    nonzero_epsilon2 = s_epsilon2(s_epsilon2 ~= 0);
    if isempty(nonzero_epsilon2)
        epsilon2_min=0;
        epsilon2_max=0;
    else
        epsilon2_min=min(nonzero_epsilon2);
        epsilon2_max=max(nonzero_epsilon2);
    end

    if i==j
        delta=sin(vector_angle)^2-a_true(j)^2-b_true(j)^2+2*abs(a_true(j))*b_true(j)*cos(vector_angle);
        delta2=sin(vector_angle)^2-a_true(j)^2-b_true(j)^2-2*abs(a_true(j))*b_true(j)*cos(vector_angle);
    else
        delta=sin(vector_angle)^2-a_true(j)^2-b_true(j)^2+2*a_true(j)*b_true(j)*cos(vector_angle);
        delta2=delta;
    end

    alpha=acos(a_true(j))*180/pi;
    beta=acos(b_true(j))*180/pi;
    alpha_beta=[abs(alpha-beta),min([alpha+beta,360-alpha-beta])];

    save_RMSE_element_meas(k,:)=sqrt(diag(P_element_meas));
    save_RMSE_element_theory(k,:)=RMSE_element;
    save_RMSE_phi_meas(k,:)=sqrt(diag(P_phi_meas));
    save_RMSE_phi_triad_meas(k,:)=sqrt(diag(P_phi_triad_meas));
    save_RMSE_phi_quest_meas(k,:)=sqrt(diag(P_phi_quest_meas));
    save_RMSE_phi_theory(k,:)=RMSE_phi;
    save_RMSE_phi_triad_theory(k,:)=RMSE_phi_triad;
    save_RMSE_phi_quest_theory(k,:)=RMSE_phi_quest;
    save_RMSE_dis_meas(k)=RMSE_dis_meas;
    save_RMSE_dis_triad_meas(k)=RMSE_dis_triad_meas;
    save_RMSE_dis_quest_meas(k)=RMSE_dis_quest_meas;
    save_RMSE_dis_theory(k)=RMSE_dis;
    save_theory(k,:)=[denom,denom2,abs(denom*denom2)];
    save_atti_true(k,:)=true_attitude;
    save_theta_true(k)=vector_angle;
    save_N_count_crit1(k,:)=N_count_crit1;
    save_N_count_crit2(k,:)=N_count_crit2;
    save_mean_dis1(k,:)=[mean_dis_atti1,mean_dis_vec1];
    save_mean_dis2(k,:)=[mean_dis_atti2,mean_dis_vec2];
    save_norm_violation_count(k)=sum(s_exc(:,1));
    save_delta_revision_count(k,:)=sum(s_exc(:,2:3),1);
    save_i_j_a_b(k,:)=[k,i,j,a_true(i),a_true(cyclicIndex(j,1)),a_true(cyclicIndex(j,2)),b_true(j),...
        delta,delta2];
    save_exc(k,:)=[sum(s_exc),max(nonzeros(s_exc_value(:,1))),min(nonzeros(s_exc_value(:,1))),...
        delta,delta2,max(nonzeros(s_exc_value(:,2))),min(nonzeros(s_exc_value(:,2))),...
        max(nonzeros(s_exc_value(:,3))),min(nonzeros(s_exc_value(:,3))),...
        epsilon_min,epsilon_max,epsilon2_min,epsilon2_max,...
        a_true(cyclicIndex(i,1)),a_true(cyclicIndex(i,2)),...
        alpha_beta,vector_angle*180/pi,min(abs(vector_angle*180/pi-alpha_beta))];
end

[~,I]=sort(save_theory(:,3),'descend');

figure(1)
s_gca=[];
set(gcf,'unit','centimeters','position',[1,1,18.1,8]);
subplot(2,1,1);
plot(save_theory(I,3),'.-','linewidth',1.4,'MarkerSize',10);
grid on
ylabel('|\tau|')
axis([1,num_test_cases,0,0.7])
set(gca,'FontName','Times New Roman','FontSize',8);
x_new = cellstr(num2str(I));
xticks(1:num_test_cases);
xticklabels(x_new);
set(gca, 'XTickLabelRotation', 45);
s_gca=[s_gca,gca];
subplot(2,1,2);
plot(save_RMSE_dis_meas(I)*180/pi,'linewidth',1.6);
hold on;
plot(save_RMSE_dis_triad_meas(I)*180/pi,'-.','linewidth',1.6);
plot(save_RMSE_dis_quest_meas(I)*180/pi,'--','linewidth',1.6);
hold off        
grid on;
xlabel('Test \it k')
ylabel('Attitude RMSE (deg)');
axis([1,num_test_cases,-inf,inf])
set(gca,'FontName','Times New Roman','FontSize',8);
legend('MIAD','TRIAD','QUEST','Location','northwest')
xticks(1:num_test_cases);
xticklabels(x_new);
set(gca, 'XTickLabelRotation', 45);
set(gca,'LooseInset', get(gca,'TightInset'))
s_gca=[s_gca,gca];
removeUnequalSubplotWhitespace(s_gca,2,[0.5,0.5])
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
print(gcf,fullfile(results_dir,'figure_1_attitude_rmse.pdf'),'-dpdf','-r0')

sorted_results = table(I,save_i_j_a_b(I,2),save_i_j_a_b(I,3), ...
    save_theory(I,3),save_RMSE_dis_meas(I)*180/pi, ...
    save_RMSE_dis_triad_meas(I)*180/pi,save_RMSE_dis_quest_meas(I)*180/pi, ...
    'VariableNames',{'TestCase','UnknownAxisI','ObservedAxisJ','AbsTau', ...
    'MiadRmseDeg','TriadRmseDeg','QuestRmseDeg'});
writetable(sorted_results,fullfile(results_dir,'simulation_results.csv'));
diagnostic_results = table((1:M)',save_i_j_a_b(:,2),save_i_j_a_b(:,3), ...
    save_i_j_a_b(:,4),save_i_j_a_b(:,5),save_i_j_a_b(:,6), ...
    save_i_j_a_b(:,7),save_i_j_a_b(:,8),save_i_j_a_b(:,9), ...
    save_N_count_crit1(:,3)/N,save_mean_dis1(:,1), ...
    save_N_count_crit2(:,3)/N,save_mean_dis2(:,1), ...
    'VariableNames',{'TestCase','UnknownAxisI','ObservedAxisJ','A_i', ...
    'A_j1','A_j2','B_j','DeltaBranchPlus','DeltaBranchMinus', ...
    'WrongRateBj1','MeanWrongErrorBj1Deg','WrongRateBj2', ...
    'MeanWrongErrorBj2Deg'});
writetable(diagnostic_results, ...
    fullfile(results_dir,'simulation_diagnostics.csv'));

norm_violation_cases = find(save_norm_violation_count > 0);
norm_violation_results = table(norm_violation_cases, ...
    save_norm_violation_count(norm_violation_cases), ...
    sum(save_exc(norm_violation_cases,16:17).^2,2), ...
    save_exc(norm_violation_cases,4),save_exc(norm_violation_cases,5), ...
    'VariableNames',{'TestCase','Count','TrueSquaredNorm', ...
    'ObservedSquaredNormMax','ObservedSquaredNormMin'});
writetable(norm_violation_results, ...
    fullfile(results_dir,'norm_violation_results.csv'));

negative_delta_test_cases = find(any(save_delta_revision_count > 0,2));
delta_test_case = zeros(2*M,1);
delta_branch = zeros(2*M,1);
delta_revision_count = zeros(2*M,1);
delta_unknown_axis_i = zeros(2*M,1);
delta_observed_axis_j = zeros(2*M,1);
true_delta = zeros(2*M,1);
observed_delta_max = zeros(2*M,1);
observed_delta_min = zeros(2*M,1);
epsilon_min = zeros(2*M,1);
epsilon_max = zeros(2*M,1);
delta_row_count = 0;
for test_index = negative_delta_test_cases'
    for branch_index = 1:2
        is_degenerate_pair = save_i_j_a_b(test_index,2) == ...
            save_i_j_a_b(test_index,3);
        has_revision = save_delta_revision_count(test_index,branch_index) > 0;
        if has_revision || is_degenerate_pair
            delta_row_count = delta_row_count + 1;
            delta_test_case(delta_row_count) = test_index;
            delta_branch(delta_row_count) = branch_index;
            delta_revision_count(delta_row_count) = ...
                save_delta_revision_count(test_index,branch_index);
            delta_unknown_axis_i(delta_row_count) = ...
                save_i_j_a_b(test_index,2);
            delta_observed_axis_j(delta_row_count) = ...
                save_i_j_a_b(test_index,3);
            true_delta(delta_row_count) = ...
                save_exc(test_index,5+branch_index);
            max_column = 6 + 2 * branch_index;
            observed_delta_max(delta_row_count) = ...
                save_exc(test_index,max_column);
            observed_delta_min(delta_row_count) = ...
                save_exc(test_index,max_column+1);
            epsilon_column = 10 + 2 * branch_index;
            epsilon_range = save_exc(test_index,epsilon_column:epsilon_column+1);
            if ~has_revision
                epsilon_range(:) = NaN;
            end
            epsilon_min(delta_row_count) = epsilon_range(1);
            epsilon_max(delta_row_count) = epsilon_range(2);
        end
    end
end
delta_test_case = delta_test_case(1:delta_row_count);
delta_branch = delta_branch(1:delta_row_count);
delta_revision_count = delta_revision_count(1:delta_row_count);
delta_unknown_axis_i = delta_unknown_axis_i(1:delta_row_count);
delta_observed_axis_j = delta_observed_axis_j(1:delta_row_count);
true_delta = true_delta(1:delta_row_count);
observed_delta_max = observed_delta_max(1:delta_row_count);
observed_delta_min = observed_delta_min(1:delta_row_count);
epsilon_min = epsilon_min(1:delta_row_count);
epsilon_max = epsilon_max(1:delta_row_count);
delta_branch_label = strings(delta_row_count,1);
delta_branch_label(delta_branch == 1) = "PositiveRoot";
delta_branch_label(delta_branch == 2) = "NegativeRoot";
negative_delta_results = table(delta_test_case,delta_revision_count, ...
    delta_unknown_axis_i,delta_observed_axis_j,delta_branch_label, ...
    true_delta,observed_delta_max,observed_delta_min, ...
    epsilon_min,epsilon_max, ...
    'VariableNames',{'TestCase','Count','UnknownAxisI','ObservedAxisJ', ...
    'DeltaBranch','TrueDelta','ObservedDeltaMax','ObservedDeltaMin', ...
    'EpsilonMin','EpsilonMax'});
writetable(negative_delta_results, ...
    fullfile(results_dir,'negative_delta_results.csv'));

category1_test_cases = [5;23;43;9;17;18;24;25;41;44;48;4];
category1_conditions = [repmat("a_i near zero",3,1); ...
    repmat("Delta near zero",8,1);"a_i and Delta near zero"];
category1_misdetection_results = table(category1_conditions, ...
    category1_test_cases,save_i_j_a_b(category1_test_cases,4), ...
    save_i_j_a_b(category1_test_cases,8), ...
    save_i_j_a_b(category1_test_cases,9), ...
    save_i_j_a_b(category1_test_cases,2), ...
    save_i_j_a_b(category1_test_cases,3), ...
    save_N_count_crit1(category1_test_cases,3)/num_monte_carlo_samples, ...
    save_mean_dis1(category1_test_cases,1), ...
    save_N_count_crit2(category1_test_cases,3)/num_monte_carlo_samples, ...
    save_mean_dis2(category1_test_cases,1), ...
    'VariableNames',{'Condition','TestCase','TrueA_i','TrueDeltaPlus', ...
    'TrueDeltaMinus', ...
    'UnknownAxisI','ObservedAxisJ','WDRUsingBj1','MeanErrorUsingBj1Deg', ...
    'WDRUsingBj2','MeanErrorUsingBj2Deg'});
writetable(category1_misdetection_results, ...
    fullfile(results_dir,'category1_misdetection_results.csv'));

category2_test_cases = [22;47;23;26;30];
category2_conditions = ["a_j1 near zero";"a_j1 near zero"; ...
    "a_j2 near zero";"b_j near zero (i = j)"; ...
    "b_j near zero (i = j)"];
category2_misdetection_results = table(category2_conditions, ...
    category2_test_cases,save_i_j_a_b(category2_test_cases,5), ...
    save_i_j_a_b(category2_test_cases,6), ...
    save_i_j_a_b(category2_test_cases,7), ...
    save_i_j_a_b(category2_test_cases,2), ...
    save_i_j_a_b(category2_test_cases,3), ...
    save_N_count_crit1(category2_test_cases,3)/num_monte_carlo_samples, ...
    save_mean_dis1(category2_test_cases,1), ...
    save_N_count_crit2(category2_test_cases,3)/num_monte_carlo_samples, ...
    save_mean_dis2(category2_test_cases,1), ...
    'VariableNames',{'Condition','TestCase','TrueA_j1','TrueA_j2', ...
    'TrueB_j','UnknownAxisI','ObservedAxisJ','WDRUsingBj1', ...
    'MeanErrorUsingBj1Deg','WDRUsingBj2','MeanErrorUsingBj2Deg'});
writetable(category2_misdetection_results, ...
    fullfile(results_dir,'category2_misdetection_results.csv'));
save(fullfile(results_dir,'simulation_results.mat'), ...
    'save_i_j_a_b','save_theory','save_RMSE_dis_meas', ...
    'save_RMSE_dis_triad_meas','save_RMSE_dis_quest_meas', ...
    'save_N_count_crit1','save_N_count_crit2','save_mean_dis1', ...
    'save_mean_dis2','save_exc','save_norm_violation_count', ...
    'save_delta_revision_count','norm_violation_results', ...
    'negative_delta_results','category1_misdetection_results', ...
    'category2_misdetection_results');
