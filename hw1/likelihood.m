function [ output_args ] = likelihood( Mu, P_x, C )
[V,D] = eig(P_x);
theta = 2*pi*[0:0.001:1];
for jj = 1:length(C),
    for ii = 1:length(theta),
        e(:,ii) = Mu + ...
            C(jj) * chol(P_x,'lower')*[cos(theta(ii));sin(theta(ii))];
    end
    axis equal, hold on
    plot(e(1,:), e(2,:), 'Color','b');
end
quiver(Mu(1),Mu(2),V(1,1),V(2,1), 'Color','k')
quiver(Mu(1),Mu(2),V(1,2),V(2,2), 'Color','k')
end

