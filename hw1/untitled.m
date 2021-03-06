axis equal, hold on


%# eigen decomposition [sorted by eigen values]
[V D] = eig( P_x);     %#' cov(X0)
[D order] = sort(diag(D), 'descend');
D = diag(D);
V = V(:, order);

t = linspace(0,2*pi,100);
e = [cos(t) ; sin(t)];        %# unit circle
VV = V*sqrt(D);               %# scale eigenvectors
e = VV*e;

%# plot cov and major/minor axes
plot(e(1,:), e(2,:), 'Color','k');
quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
