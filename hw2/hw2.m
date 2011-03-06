%% Problem 3
% Use least squares to identify a gyroscopes (a) and (b).
a = 0.1; b = 0.5;
sigma = 0.3;

omega = 2*pi;
sim_length = 1000;
t = linspace(0,1,sim_length);

for ii = 1:1000
    r = 100 * sin(omega * t);
    g = a * r + b + sigma*randn(sim_length,1)';
    
    H = [r', ones(length(r),1)];
    x(ii,:) = (pinv(H) * g')';
end

x_mean = mean(x)
x_std = sqrt(var(x))

%%

a = 0.1; b = 0.5;
sigma = 0.3;

omega = 2*pi;
t = linspace(0.1,1,1000);

r = 100 * sin(omega * t);
g = a * r + b + sigma*randn(1000,1)';

P = 10000 * eye(2);
x(:,1) = 0;

for ii=2:length(r),
    pi_mat = r(ii) * P;
    gamma = pi_mat * r(ii);
    k(ii) = pi_mat'/gamma;
    alpha(ii) = g(ii) - x(:,ii-1)'*g(ii);
    x(:,ii) = x(:,ii-1) + k(ii)*alpha(ii);
    Pprime = k(ii)*pi_mat;
    P = P - Pprime;
end

%% Problem 4

x = [];

for ii=1:10
    % Simulate the discrete transfer function.
    numd = 0.25 * [1 -0.8];
    dend = [1 -1.9 0.95];
    u = randn(1000,1);
    y = dlsim(numd,dend,u);
    sigma = 1;
    Y = y + sigma*randn(1000,1);

    % Form the H-matrix
    H = [u(2:end-1), u(1:end-2), -Y(2:end-1), -Y(3:end)];

    % Determine x hat.
    x(ii,:) = (pinv(H) * Y(1:end-2))';
    
    y_hat = H*x(ii,:)';
    
    r = Y(1:end-2) - y_hat;
end

mean(x)
sqrt(var(x))

numd_lsqr = [x(ii,1) x(ii,2)];
dend_lsqr = [1 x(ii,3) x(ii,4)];

orig = tf(numd,dend,1/1000);
lsqrs = tf(numd_lsqr,dend_lsqr,1/1000);

bode(lsqrs)