%%
t = [0:0.001:0.5];

x = zeros(2,length(t));
y = zeros(1,length(t));

x_hat = zeros(2,length(t));
y_hat = zeros(1,length(t));

u = ones(1,length(t));
r = zeros(1,length(t));

for(n = [2:length(t)-1]),
   
   r(n) = u(n) - y_hat(n-1);
   
   y(n) = sys_d.C * x(:,n);
   x(:,n+1) = sys_d.A * x(:,n) + sys_d.B * r(n);
   
   y_hat(n) = -K_d * x_hat(:,n);
   x_hat(:,n+1) = Ahat_d * x_hat(:,n) + ...
       sys_d.B * r(n) + ...
       L_d * (y(n)-y_hat(n));
   
   
   
end