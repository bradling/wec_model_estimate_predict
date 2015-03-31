function [varargout] = plot_rad_approx(obj)

h = figure;
set(h, 'color', 'w')
set(h, 'units', 'inches')
set(h, 'position', [1 1 6.5 6])
set(h, 'name', 'S.S. Rad. Approx. Impulse Responses')


tlims = obj.radIrf.t([1 end]);
titleStrings = {'Surge', 'Heave', 'Pitch', 'Surge-Heave', ...
                'Surge-Pitch', 'Heave-Pitch'};
            
for ii = 1:6
    subplot(3, 2, ii)
    [y, t] = impulse(obj.ssRad{ii});
    plot(obj.radIrf.t, obj.radIrf.irf(ii,:), t, y, 'r--')
    set(gca, 'fontsize', 9)
    xlim(tlims)
    grid on
    xlabel('Time (sec)')
    ylabel('Rad. Impulse Resp.')
    legend('True', 'Approx.', 'location', 'best')
    title(titleStrings{ii})
end

if nargout > 0
    varargout{1} = h;
end