function vol = body_volume(s)
% body_volume gives the 3D volume of the axis-symmetric body
% s - quadrature of the generating arc (from top to bottom, half-circle)

% add a negative sign because the curve is going downward, so dy/dt flips
% Vy = int_{a}^{b} pi x^2 (dy/dt) dt 
% https://en.wikipedia.org/wiki/Solid_of_revolution

vol = -pi*s.w'*((real(s.x).^2).*imag(s.xp));
end