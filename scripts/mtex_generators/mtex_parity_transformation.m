function computed = mtex_parity_transformation(currentCase)
%MTEX_PARITY_TRANSFORMATION Orientation-relationship parity outputs from MTEX.
%
%   Handles two operations of the `or_transformation_v1` campaign:
%
%     or_misorientation_representative
%       The relationship as a symmetry-reduced misorientation (angle, axis)
%       plus its variant count.
%
%     or_fit_from_orientation_pairs
%       Child orientations are DERIVED from the listed parents and variant
%       indices of a generating relationship, then calcParent2Child is seeded
%       from a nominal relationship and must recover the generating one.
%
%   IMPORTANT: this handler has NOT been executed against a real MTEX
%   installation — none was available on the machine where the campaign was
%   authored. The PyTex side of the campaign is generated and verified; this
%   file encodes the intended MTEX comparison and should be expected to need
%   small API corrections on first run. Do not treat a failure here as a
%   PyTex/MTEX disagreement until the script itself runs cleanly. Record any
%   corrections in scripts/mtex_generators/README.md.

[csParent, SS] = mtex_phase_symmetry(currentCase.phase);
input = currentCase.input;
csChild = mtex_phase_symmetry(input.child_phase);

switch string(currentCase.operation)
    case "or_misorientation_representative"
        p2c = local_named_relationship(string(input.relationship), csParent, csChild);
        [angleDeg, axisSorted] = local_misorientation_descriptor(p2c);
        computed = struct();
        computed.results = struct( ...
            "relationship", string(input.relationship), ...
            "misorientation_angle_deg", angleDeg, ...
            "misorientation_axis_sorted_abs", axisSorted, ...
            "misorientation_quaternion_wxyz", mtex_quaternion_wxyz(p2c), ...
            "variant_count", length(variants(p2c)));
        computed.artifacts = [];
        computed.notes = ["MTEX parent-to-child misorientation and variant count."];

    case "or_fit_from_orientation_pairs"
        nominal = local_named_relationship(string(input.nominal_relationship), csParent, csChild);
        generating = local_named_relationship(string(input.generating_relationship), csParent, csChild);

        euler = input.parent_euler_deg;
        parentOri = orientation.byEuler( ...
            euler(:, 1) * degree, euler(:, 2) * degree, euler(:, 3) * degree, csParent, SS);

        % Child orientations through the generating relationship. MTEX's
        % variants(p2c) returns the child-to-parent variant set; a child is
        % obtained by composing the parent with the inverse variant rotation,
        % matching the canonical C = P V' composition used by PyTex.
        generatingVariants = variants(generating);
        indices = input.generating_variant_indices(:);
        childOri = parentOri;
        for k = 1:numel(indices)
            childOri(k) = parentOri(k) .* inv(generatingVariants(indices(k))); %#ok<MINV>
        end

        % calcParent2Child refines the parent-to-child rotation from measured
        % parent/child pairs, seeded from the nominal relationship.
        measured = inv(parentOri) .* childOri;
        fitted = calcParent2Child(measured, nominal);

        [angleDeg, axisSorted] = local_misorientation_descriptor(fitted);
        residuals = angle(measured, fitted) / degree;

        computed = struct();
        computed.results = struct( ...
            "nominal_relationship", string(input.nominal_relationship), ...
            "generating_relationship", string(input.generating_relationship), ...
            "pair_count", numel(indices), ...
            "child_euler_deg", mtex_euler_bunge_deg(childOri), ...
            "fitted_angle_deg", angleDeg, ...
            "fitted_axis_sorted_abs", axisSorted, ...
            "fitted_quaternion_wxyz", mtex_quaternion_wxyz(fitted), ...
            "deviation_from_nominal_deg", angle(fitted, nominal) / degree, ...
            "mean_residual_deg", mean(residuals(:)), ...
            "max_residual_deg", max(residuals(:)));
        computed.artifacts = [];
        computed.notes = [ ...
            "MTEX calcParent2Child fit seeded from the nominal relationship."; ...
            "Child orientations are derived from MTEX's own variant enumeration, so a" + ...
            " differing variant order from PyTex changes per-pair residual assignment" + ...
            " but must not change the fitted rotation."];

    otherwise
        error("PyTexParity:UnknownTransformationOperation", ...
            "Unsupported transformation operation: %s", currentCase.operation);
end
end


function p2c = local_named_relationship(name, csParent, csChild)
%LOCAL_NAMED_RELATIONSHIP Map the campaign's relationship name onto MTEX.
switch name
    case "bain"
        p2c = orientation.Bain(csParent, csChild);
    case "burgers"
        p2c = orientation.Burgers(csParent, csChild);
    case "greninger_troiano"
        p2c = orientation.GreningerTrojano(csParent, csChild);
    case "kurdjumov_sachs"
        p2c = orientation.KurdjumovSachs(csParent, csChild);
    case "nishiyama_wassermann"
        p2c = orientation.NishiyamaWassermann(csParent, csChild);
    case "pitsch"
        p2c = orientation.Pitsch(csParent, csChild);
    otherwise
        error("PyTexParity:UnknownRelationship", "Unknown relationship: %s", name);
end
end


function [angleDeg, axisSorted] = local_misorientation_descriptor(p2c)
%LOCAL_MISORIENTATION_DESCRIPTOR Angle and point-group-invariant axis descriptor.
%
%   The symmetry-reduced representative is defined only up to the crystal
%   point group, so the axis is reported as its absolute components sorted
%   descending — the comparison convention declared by the campaign.
angleDeg = angle(p2c) / degree;
axisVector = axis(p2c);
components = [axisVector.x, axisVector.y, axisVector.z];
components = components ./ norm(components);
axisSorted = sort(abs(components), "descend");
end
