function computed = mtex_parity_miller(currentCase)
%MTEX_PARITY_MILLER Compute Miller geometry and symmetry outputs with MTEX.

[CS, ~] = mtex_phase_symmetry(currentCase.phase);
input = currentCase.input;

if isfield(input, "plane_hkil_a")
    pA = Miller(input.plane_hkil_a(1), input.plane_hkil_a(2), input.plane_hkil_a(3), input.plane_hkil_a(4), CS);
    pB = Miller(input.plane_hkil_b(1), input.plane_hkil_b(2), input.plane_hkil_b(3), input.plane_hkil_b(4), CS);
    dA = Miller(input.direction_uvtw_a(1), input.direction_uvtw_a(2), input.direction_uvtw_a(3), input.direction_uvtw_a(4), CS, "uvw");
    dB = Miller(input.direction_uvtw_b(1), input.direction_uvtw_b(2), input.direction_uvtw_b(3), input.direction_uvtw_b(4), CS, "uvw");
else
    pA = Miller(input.plane_hkl_a(1), input.plane_hkl_a(2), input.plane_hkl_a(3), CS);
    pB = Miller(input.plane_hkl_b(1), input.plane_hkl_b(2), input.plane_hkl_b(3), CS);
    dA = Miller(input.direction_uvw_a(1), input.direction_uvw_a(2), input.direction_uvw_a(3), CS, "uvw");
    dB = Miller(input.direction_uvw_b(1), input.direction_uvw_b(2), input.direction_uvw_b(3), CS, "uvw");
end

results = struct();
results.plane_plane_angle_deg = double(angle(pA, pB) / degree);
results.direction_direction_angle_deg = double(angle(dA, dB) / degree);
results.plane_a_family_hkl = mtex_miller_array_indices(symmetrise(pA));
results.direction_a_family_uvw = mtex_miller_array_indices(symmetrise(dA));

try
    results.plane_a_d_spacing_angstrom = double(norm(pA));
catch
    results.plane_a_d_spacing_angstrom = NaN;
end

if isfield(input, "plane_hkil_a")
    results.plane_a_hkl = mtex_miller_indices(pA);
    results.direction_a_uvw = mtex_miller_indices(dA);
end

computed = struct();
computed.results = results;
computed.artifacts = [];
computed.notes = ["MTEX Miller geometry result."];
end

