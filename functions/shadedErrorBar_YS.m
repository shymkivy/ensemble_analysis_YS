function s = shadedErrorBar_YS(plot_time, temp_trace_mean,temp_trace_sem, col)


s = shadedErrorBar(plot_time, temp_trace_mean,temp_trace_sem);
if exist('col', 'var')
    s.patch.FaceColor = col;
    s.mainLine.Color = col;
    s.edge(1).Color = min(col + 0.5,1);
    s.edge(2).Color = min(col + 0.5,1);
end



end