%% MECH7710 - HW1
% Random Variables and Probability
%% Problem 1

%%% Part A - 6 dice numbered 1,2,3,4,5,6
pdf_matrix = 1/6 * ones(6,6);
pdf_combined = nfoldconv(pdf_matrix);
[mean,variance] = statistics(pdf_combined);
true_norm = normpdf([1:length(pdf_combined)],mean,sqrt(variance));
if plots,
    createfigure(true_norm, pdf_combined, 'Problem 1 - Part A')
end

%%% Part B - 6 dice numbered 4,5,6,7,8,9
pdf_matrix = 1/6 * [zeros(6,3),ones(6,6)];
pdf_combined = nfoldconv(pdf_matrix);
[mean,variance] = statistics(pdf_combined);
true_norm = normpdf([1:length(pdf_combined)],mean,sqrt(variance));
if plots,
    createfigure(true_norm, pdf_combined, 'Problem 1 - Part B')
end

%%% Part C - 6 dice numbered 1,1,3,3,3,5
pdf_matrix = [...
    1/3 * ones(6,1), zeros(6,1), ...
    1/2 * ones(6,1), zeros(6,1), ...
    1/6 * ones(6,1), zeros(6,1)];
pdf_combined = nfoldconv(pdf_matrix);
[mean,variance] = statistics(pdf_combined);
true_norm = normpdf([1:length(pdf_combined)],mean,sqrt(variance));
if plots,
    createfigure(true_norm, pdf_combined, 'Problem 1 - Part C')
end

%%% Part D - 3 dice numbered 1,2,3,4,5,6 and 3 numbered 1,1,3,3,3,5
pdf_matrix = [...
    1/6*ones(3,6); ...
    1/3 * ones(3,1), zeros(3,1), ...
    1/2 * ones(3,1), zeros(3,1), ...
    1/6 * ones(3,1), 1/6*zeros(3,1)];

pdf_combined = nfoldconv(pdf_matrix);
[mean,variance] = statistics(pdf_combined);
true_norm = normpdf([1:length(pdf_combined)],mean,sqrt(variance));
if plots,
    createfigure(true_norm, pdf_combined, 'Problem 1 - Part D')
end

%% Problem 2

%%% Part A - Mean, Central Moment, Mean Squared, Variance and Covariance
pdf_die_1 = 1/6 * ones(1,6);
pdf_die_2 = 1/6 * ones(1,6);

die = 1:6;

joint_pdf = pdf_die_1' * pdf_die_2;

[mean_1, variance_1, c_moment_1, mean_sq_1] = statistics(pdf_die_1);

% Covariance of 2 independant variables is 0.
covariance = sum(sum((die-mean_1)'*(die-mean_1) * joint_pdf)); 


%%% Part B - Covariance Matrix
P = [variance_1, covariance; covariance, variance_1];

%% Problem 6
P_x = [2, 1; 1, 4];

[V,D] = eig(P_x);
[D order] = sort(diag(D), 'descend');
D = diag(D);
V = V(:, order);

theta = 2*pi*[0:0.001:1];
for ii = 1:length(theta),
    e(:,ii) = chol(P_x,'lower')*[cos(theta(ii));sin(theta(ii))];
end

axis equal, hold on

C = [0.25,1,1.25];

for ii = 1:length(C),
    plot(C(ii)*e(1,:), C(ii)*e(2,:), 'Color','b');
end

quiver(0,0,V(1,1),V(2,1), 'Color','k')
quiver(0,0,V(1,2),V(2,2), 'Color','k')

