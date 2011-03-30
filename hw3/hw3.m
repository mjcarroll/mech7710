%% Problem 1

x = [];

for ii=1:10
    % Simulate the discrete transfer function.
    numd = 0.25 * [1 -0.8];
    dend = [1 -1.9 0.95];
    u = randn(1000,1);
    y = dlsim(numd,dend,u);
    sigma = 0.01;
    Y = y + sigma*randn(1000,1);

    % Form the H-matrix
    H = [u(2:end-1), u(1:end-2), -Y(2:end-1), -Y(3:end)];

    % Determine x hat.
    x(ii,:) = (pinv(H) * Y(1:end-2))';
    y_hat = H*x(ii,:)';
    r = Y(1:end-2) - y_hat;
end

mean(x);
sqrt(var(x));

numd_lsqr = [x(ii,1) x(ii,2)];
dend_lsqr = [1 x(ii,3) x(ii,4)];

orig = tf(numd,dend,1/1000);
lsqrs = tf(numd_lsqr,dend_lsqr,1/1000);

bode(lsqrs)
hold on
bode(orig)

