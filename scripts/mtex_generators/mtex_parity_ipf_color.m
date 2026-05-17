function computed = mtex_parity_ipf_color(currentCase)
%MTEX_PARITY_IPF_COLOR Compute MTEX IPF RGB colors from Euler orientations.

[CS, SS] = mtex_phase_symmetry(currentCase.phase);
input = currentCase.input;
angles = input.euler_deg;
ori = orientation.byEuler(angles(:, 1) * degree, angles(:, 2) * degree, angles(:, 3) * degree, CS, SS);

key = ipfHSVKey(CS);
key.inversePoleFigureDirection = mtex_specimen_direction(input.specimen_direction);
rgb = key.orientation2color(ori);

crystalDirections = inv(ori) * key.inversePoleFigureDirection;
reduced = project2FundamentalRegion(crystalDirections, CS);

computed = struct();
computed.results = struct( ...
    "rgb", rgb, ...
    "crystal_directions", mtex_vector_array_xyz(crystalDirections), ...
    "sector_reduced_directions", mtex_vector_array_xyz(reduced), ...
    "color_space", "srgb_0_1", ...
    "specimen_direction", mtex_vector_xyz(key.inversePoleFigureDirection));
computed.artifacts = [];
computed.notes = ["MTEX IPF color result. RGB tolerances should account for color-key implementation details."];
end

