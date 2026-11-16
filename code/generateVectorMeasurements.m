function [observed_a_samples,observed_b_samples,true_a,true_b,true_attitude] = ...
    generateVectorMeasurements(unknown_index_i,observed_index_j,sigma_a_i1, ...
    sigma_a_i2,sigma_b_j,sigma_a_i,sigma_b_j1,sigma_b_j2,reference_p, ...
    reference_q,num_samples)
%GENERATEVECTORMEASUREMENTS Generate noisy sensor-frame vector components.

true_attitude = [180 * (2 * rand - 1); ...
    90 * (2 * rand - 1); 180 + 180 * (2 * rand - 1)] * pi / 180;
Cnb = eulerZxyToCnb(true_attitude);
Cbn = Cnb';
true_a = Cbn * reference_p;
true_b = Cbn * reference_q;
observed_a_samples = zeros(num_samples,3);
observed_b_samples = zeros(num_samples,3);

for sample_index = 1:num_samples
    observed_a = zeros(1,3);
    observed_a(cyclicIndex(unknown_index_i,1)) = ...
        true_a(cyclicIndex(unknown_index_i,1)) + sigma_a_i1 * randn;
    observed_a(cyclicIndex(unknown_index_i,2)) = ...
        true_a(cyclicIndex(unknown_index_i,2)) + sigma_a_i2 * randn;
    observed_a(unknown_index_i) = ...
        true_a(unknown_index_i) + sigma_a_i * randn;

    observed_b = zeros(1,3);
    observed_b(observed_index_j) = ...
        true_b(observed_index_j) + sigma_b_j * randn;
    observed_b(cyclicIndex(observed_index_j,1)) = ...
        true_b(cyclicIndex(observed_index_j,1)) + sigma_b_j1 * randn;
    observed_b(cyclicIndex(observed_index_j,2)) = ...
        true_b(cyclicIndex(observed_index_j,2)) + sigma_b_j2 * randn;

    observed_a_samples(sample_index,:) = observed_a;
    observed_b_samples(sample_index,:) = observed_b;
end
end
