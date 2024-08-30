function val_target = interpolate_by_panel(s, val, t_target)
% interpolate for quantities on the boundary
% val is an attribute of s, for example, val=s.xp
% t_target is any give vector (monotonic increasing) in the domain for t

val_target = zeros(size(t_target));

for pnum = 1:numel(s.tlo)
    tlo = s.tlo(pnum); thi = s.thi(pnum);
    % find t_inpanel
    t_ind = find(t_target>=tlo&t_target<=thi);
    if ~isempty(t_ind)
        t_inpanel = t_target(t_ind);
        t_inpanel_rescale = - 1 + 2*(t_inpanel-tlo)/(thi-tlo);
        val_target(t_ind) = interpmat_1d(t_inpanel_rescale,gauss(s.p)) * val(((pnum-1)*s.p+1):(pnum*s.p));
    end    
end

end