function createfigure(true_norm1, pdf_combined1, figure_name)

% Create figure
figure1 = figure('XVisual',...
    '0x27 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)',...
    'Name',figure_name);
% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');
% Create plot
plot(true_norm1,'Parent',axes1,'DisplayName','Normal Distribution');
% Create stem
stem(pdf_combined1,'Parent',axes1,'DisplayName','Combined PDF Function');
% Create xlabel
xlabel('Sum of 6 Dice');
% Create ylabel
ylabel('Probability Distribution');
% Create title
title('Probability Distribution Function of Dice Roll');
% Create legend
legend(axes1,'show');

end