subnames = {};
for i = 31522:31767, subnames{i-31521}=sprintf('%0.6i', i); end
sub_outi={};
idx_out.all=ismember(cellfun(@(x) x(1:6),subnames,'UniformOutput',false),sub_outi);
subnames=(subnames(~idx_out.all));

