function computed = mtex_parity_odf(currentCase)
%MTEX_PARITY_ODF Compute sampled ODF parity outputs with MTEX.

if strcmp(string(currentCase.operation), "odf_from_xrdml")
    error("PyTexParity:PendingXRDML", "XRDML ODF cases are pending until fixture files are provided.");
end

[CS, SS] = mtex_phase_symmetry(currentCase.phase);
input = currentCase.input;
support = input.support_euler_deg;
query = input.query_euler_deg;
weights = input.weights(:);
normalizedWeights = weights ./ sum(weights);

supportOri = orientation.byEuler(support(:, 1) * degree, support(:, 2) * degree, support(:, 3) * degree, CS, SS);
queryOri = orientation.byEuler(query(:, 1) * degree, query(:, 2) * degree, query(:, 3) * degree, CS, SS);

halfwidth = input.kernel.halfwidth_deg * degree;
kernel = deLaValleePoussinKernel("halfwidth", halfwidth);
odf = calcDensity(supportOri, "weights", weights, "kernel", kernel);
density = eval(odf, queryOri);

computed = struct();
computed.results = struct( ...
    "odf_representation", "sampled_density", ...
    "support_euler_deg", support, ...
    "query_euler_deg", query, ...
    "normalized_weights", normalizedWeights, ...
    "query_density", density(:), ...
    "summary", struct( ...
        "min", min(density(:)), ...
        "max", max(density(:)), ...
        "mean", mean(density(:))));
computed.artifacts = [];
computed.notes = ["MTEX sampled ODF result from weighted Euler orientations."];
end

