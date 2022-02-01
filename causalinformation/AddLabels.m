function AddLabels(D0,clrs)
figure;
hold on
% Plot whatever you like
x = 1:10;
y = NaN;
for i=1:length(D0)
    plot(x, y .* x,clrs{i},'LineWidth',2)
end

leg = string(D0);
for l = 1:length(D0)
   leg(l) = 'D(t) = ' + leg(l);
end
% Initial values to capture the entire legend
% Should fit most modern screens
set(gcf,'Position',[0,0,1024,1024]);
% Call the legend to your choice, I used a horizontal legend here
legend_handle = legend(leg,'Orientation','horizontal');
% Set the figure Position using the normalized legend Position vector
% as a multiplier to the figure's current position in pixels
% This sets the figure to have the same size as the legend
set(gcf,'Position',(get(legend_handle,'Position')...
    .*[0, 0, 1, 1].*get(gcf,'Position')));
% The legend is still offset so set its normalized position vector to
% fill the figure
set(legend_handle,'Position',[0,0,1,1]);
% Put the figure back in the middle screen area
set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);
end