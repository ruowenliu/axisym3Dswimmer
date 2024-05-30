function sa = surface_area_prolate(a,b)
% calculate surface area for prolate
% a = semi axis for x-axis
% b = semi axis for z-axis

if a==b
    sa = 4*pi*a*b;
else
    ec = sqrt(1-a^2/b^2); % eccentricity
    sa = 2*pi*a^2*(1+b/a/ec*asin(ec));
end
end
