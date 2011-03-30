a = 0.1; b = 0.5;
sigma = 0.3;

omega = 2*pi;
sim_length = 10;
t = linspace(0,1,sim_length);

r = 100 * sin(omega * t);
g = a * r + b + sigma*randn(sim_length,1)';

H = [r', ones(length(r),1)];
x = (pinv(H) * g')';
