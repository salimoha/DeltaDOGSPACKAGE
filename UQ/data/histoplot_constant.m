function [h]=histoplot_constant(y,N,lb,ub,index)

intervalWidth = (ub - lb)/N;
x = lb:intervalWidth:ub;

ncount = histc(y,x);
%relativefreq = ncount/length(y);
figure(index)
h=bar(x-intervalWidth/2, ncount,1);
xlim([lb ub])
set(gca,   'FontSize', 18);
%set(gca, 'xtick', x)
end