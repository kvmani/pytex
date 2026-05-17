function angles = mtex_euler_bunge_deg(orientationLike)
%MTEX_EULER_BUNGE_DEG Return Bunge Euler angles in degrees.
[phi1, Phi, phi2] = Euler(orientationLike);
angles = double([phi1, Phi, phi2] ./ degree);
end

