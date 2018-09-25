function plot_base_graph(BG)
    V = get_3gpp_base_graph(BG,0);
    figure
    colormap([1 1 1; 0 0 0]);
    image(0:size(V,2)-1,0:size(V,1)-1,(V>0));
    title(sprintf('Base graph %d',BG));
    xlabel('column');
    ylabel('row');