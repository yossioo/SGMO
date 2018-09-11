function p_out = project_point_on_line(p,line)

pA = line(1,:);
v = diff(line) ;
v = v./norm(v);
l = v* (p-pA)';
% l = dot(v, p-pA);
p_out = pA + l(:).*repmat(v,numel(l),1);

end

