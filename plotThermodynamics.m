function plotThermodynamics(y, xlabels, name)
f = figure('visible','off');
    plot(0:numel(y), cumsum([0 y]), '-o', 'MarkerFaceColor', 'b');
    title('dG')
    ax = gca;
    ax.XTick = 0:numel(y) ;
    xlabel('KEGG reaction ID');
    set(ax,'XTickLabel',[' '; xlabels]);
    
    ylabel('dG');
    
   saveas(f,name,'pdf')
   clear f;
   clear ax;
end