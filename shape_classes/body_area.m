function val = body_area(s)
% body_area gives the 3D surface area of the axis-symmetric body
% s - quadrature of the generating arc (from top to bottom, half-circle)
val = 2*pi*s.ws'*real(s.x);
end