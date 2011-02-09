%% MECH7710 - HW1
% Random Variables and Probability
%% Problem 1

%%% Part A - 6 dice numbered 1,2,3,4,5,6
pdf_matrix = 1/6 * ones(6,6);
pdf_combined = nfoldconv(pdf_matrix);
[mean,variance] = statistics(pdf_combined)
true_norm = normpdf([1:length(pdf_combined)],mean,sqrt(variance));
if plots,
    createfigure(true_norm, pdf_combined, 'Problem 1 - Part A')
end

%%% Part B - 6 dice numbered 4,5,6,7,8,9
pdf_matrix = 1/6 * [zeros(6,3),ones(6,6)];
pdf_combined = nfoldconv(pdf_matrix);
[mean,variance] = statistics(pdf_combined)
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
[mean,variance] = statistics(pdf_combined)
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
[mean,variance] = statistics(pdf_combined)
true_norm = normpdf([1:length(pdf_combined)],mean,sqrt(variance));
if plots,
    createfigure(true_norm, pdf_combined, 'Problem 1 - Part D')
end

clear mean;
%% Problem 2

%%% Part A - Mean, Central Moment, Mean Squared, Variance and Covariance
pdf_die_1 = 1/6 * ones(1,6);
pdf_die_2 = 1/6 * ones(1,6);

die = 1:6;
joint_pdf = pdf_die_1' * pdf_die_2;

[mean_1, variance_1, c_moment_1, mean_sq_1] = statistics(pdf_die_1);

% Covariance of 2 independant variables is 0.
covariance = sum(sum((die-mean_1)'*(die-mean_1) * joint_pdf))

%%% Part B - Covariance Matrix
P = [variance_1, covariance; covariance, variance_1]

%%% Part C - Find the PDF matrix for the variables $v_1=x_1$ 
% and $v_2 = x_1 + x_2$. Use [1:6, 1:12], with the first column zeros.

v_1 = [1:6]; v_2 = [1:12];

joint_pdf = zeros(6,12); joint_v = zeros(6,12);

for ii=1:6,
    joint_pdf(ii,ii+[1:6]) = 1/36;
    for jj=1:6,
        joint_v(ii,ii+jj) = ii+jj;
    end
end

[mean_v1, variance_v1, c_moment_v1, mean_sq_v1] = statistics(pdf_die_1)

mean_v2 = sum(sum(joint_v .* joint_pdf))
mean_sq_v2 = sum(sum(joint_v.^2 .* joint_pdf))
c_moment_v2 = sum(sum((joint_v - mean_v2) .* joint_pdf))
variance_v2 = sum(sum((joint_v - mean_v2).^2 .* joint_pdf))

covariance_12 = sum(sum((v_1 - mean_v1)'*(v_2 - mean_v2) .* joint_pdf))

P = [variance_v1, covariance_12;
    covariance_12, variance_v2]
rho_12 = covariance_12/(sqrt(P(1,1)) * sqrt(P(2,2)))

%% Problem 3
% Worked on paper 1/27Th

%% Problem 4

v_0 = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5];
pdf_0 = 1/6 * ones(1,6);
pdf_new = conv(pdf_0,pdf_0);
v_new = [-5:5];
mean_v0 = sum(pdf_new .* v_new)
variance_v0 = sum((v_new - mean_v0).^2 .* pdf_new)

%%

r = 0.1; v = 0;
var_v0 = r*variance_v0

    for ii=1:1000,
        v(ii+1) = (1-r)*v(ii) + r*(v_0(randi(6)) + v_0(randi(6)));
    end
    mean_vn = mean(v)
    var_vn = var(v)

R = conv(v,fliplr(v))/(length(v)+1);
ndt = (-(length(v)-1):(length(v)-1))* 1/1000;

R_est =r .* exp(-(1-r).*abs(ndt)) .* sqrt(variance_v0);

%P_x = cov(v(1000:2000),v(1000:2000))

%scatter(v(1000:2000),v(100:1100))
%likelihood([0;0],P_x,[0.25,1.0,2.0])

%% Problem 5
% Worked on paper 1/27Th

%% Problem 6
P_x = [2, 1; 1, 4];

Mu = [0;0];
[V,D] = eig(P_x);

likelihood(Mu,P_x,[0.25, 1.0, 1.5]);

%% Problem 7
% To work on Paper

X = linspace(0,5,1000);

sigma = 2;

fx_pdf = normpdf(linspace(-5,5,1000),0,sigma);
fy_pdf = 1./(sigma * sqrt(4*pi*X)) .* exp(-X./(4*sigma^2));

figure7 = figure('Name','Problem 7');
if plots,
plot(X,fy_pdf)
hold
plot(linspace(-5,5,1000),fx_pdf)
end



