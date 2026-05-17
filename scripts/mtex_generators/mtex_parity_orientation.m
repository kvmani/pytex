function computed = mtex_parity_orientation(currentCase, mode)
%MTEX_PARITY_ORIENTATION Compute orientation parity outputs with MTEX.

[CS, SS] = mtex_phase_symmetry(currentCase.phase);
input = currentCase.input;

switch string(mode)
    case "euler"
        ori = orientation.byEuler(input.euler_deg(1) * degree, input.euler_deg(2) * degree, input.euler_deg(3) * degree, CS, SS);
    case "axis_angle"
        rot = rotation.byAxisAngle(vector3d(input.axis(1), input.axis(2), input.axis(3)), input.angle_deg * degree);
        ori = orientation(rot, CS, SS);
    case "quaternion"
        q = input.quaternion_wxyz;
        rot = rotation(quaternion(q(1), q(2), q(3), q(4)));
        ori = orientation(rot, CS, SS);
    case "matrix"
        rot = rotation.byMatrix(input.rotation_matrix);
        ori = orientation(rot, CS, SS);
    case "miller"
        plane = Miller(input.plane_hkl(1), input.plane_hkl(2), input.plane_hkl(3), CS);
        direction = Miller(input.direction_uvw(1), input.direction_uvw(2), input.direction_uvw(3), CS, "uvw");
        ori = orientation("Miller", plane, direction, CS, SS);
    case "operations"
        left = orientation.byEuler(input.left_euler_deg(1) * degree, input.left_euler_deg(2) * degree, input.left_euler_deg(3) * degree, CS, SS);
        right = orientation.byEuler(input.right_euler_deg(1) * degree, input.right_euler_deg(2) * degree, input.right_euler_deg(3) * degree, CS, SS);
        composed = left * right;
        inverseLeft = inv(left);
        computed = struct();
        computed.results = struct( ...
            "left_quaternion_wxyz", mtex_quaternion_wxyz(left), ...
            "right_quaternion_wxyz", mtex_quaternion_wxyz(right), ...
            "composed_quaternion_wxyz", mtex_quaternion_wxyz(composed), ...
            "inverse_left_quaternion_wxyz", mtex_quaternion_wxyz(inverseLeft), ...
            "misorientation_angle_deg", double(angle(left, right) / degree), ...
            "mapped_test_vector", mtex_vector_xyz(left * vector3d(input.test_vector_crystal(1), input.test_vector_crystal(2), input.test_vector_crystal(3))));
        computed.artifacts = [];
        computed.notes = ["MTEX orientation operation result."];
        return
    otherwise
        error("PyTexParity:UnknownOrientationMode", "Unknown orientation mode: %s", mode);
end

results = struct( ...
    "quaternion_wxyz", mtex_quaternion_wxyz(ori), ...
    "rotation_matrix", mtex_rotation_matrix(ori), ...
    "euler_bunge_deg", mtex_euler_bunge_deg(ori), ...
    "axis_angle", mtex_axis_angle_struct(ori), ...
    "symmetry_equivalent_count", numel(CS.rot));

if isfield(input, "map_crystal_vectors")
    mapped = zeros(size(input.map_crystal_vectors));
    for idx = 1:size(input.map_crystal_vectors, 1)
        vec = input.map_crystal_vectors(idx, :);
        mapped(idx, :) = mtex_vector_xyz(ori * vector3d(vec(1), vec(2), vec(3)));
    end
    results.mapped_crystal_vectors = mapped;
end

if strcmp(string(mode), "miller")
    results.mapped_plane_normal = mtex_vector_xyz(ori * plane);
    results.mapped_direction = mtex_vector_xyz(ori * direction);
end

computed = struct();
computed.results = results;
computed.artifacts = [];
computed.notes = ["MTEX orientation construction result."];
end

