function result = mtex_run_single_parity_case(campaign, currentCase, caseFile)
%MTEX_RUN_SINGLE_PARITY_CASE Dispatch one shared PyTex/MTEX parity case.

if ~strcmp(currentCase.status, "active")
    result = mtex_base_result(campaign, currentCase, caseFile, "skipped");
    result.results = struct("reason_pending", get_optional_string(currentCase, "reason_pending"));
    result.notes = ["Case is pending and was not computed by MTEX."];
    return
end

switch string(currentCase.operation)
    case "orientation_from_euler"
        computed = mtex_parity_orientation(currentCase, "euler");
    case "orientation_from_axis_angle"
        computed = mtex_parity_orientation(currentCase, "axis_angle");
    case "orientation_from_quaternion"
        computed = mtex_parity_orientation(currentCase, "quaternion");
    case "orientation_from_matrix"
        computed = mtex_parity_orientation(currentCase, "matrix");
    case "orientation_from_miller"
        computed = mtex_parity_orientation(currentCase, "miller");
    case "orientation_operations"
        computed = mtex_parity_orientation(currentCase, "operations");
    case "ipf_colors_from_euler"
        computed = mtex_parity_ipf_color(currentCase);
    case "miller_metrics"
        computed = mtex_parity_miller(currentCase);
    case "odf_from_discrete_euler"
        computed = mtex_parity_odf(currentCase);
    case "pole_figure_from_xrdml"
        computed = mtex_parity_pole_figure(currentCase);
    case "odf_from_xrdml"
        computed = mtex_parity_odf(currentCase);
    case "or_misorientation_representative"
        computed = mtex_parity_transformation(currentCase);
    case "or_fit_from_orientation_pairs"
        computed = mtex_parity_transformation(currentCase);
    otherwise
        error("PyTexParity:UnknownOperation", "Unsupported parity operation: %s", currentCase.operation);
end

result = mtex_base_result(campaign, currentCase, caseFile, "active");
result.results = computed.results;
result.artifacts = computed.artifacts;
result.notes = computed.notes;
end

function result = mtex_failed_result(campaign, currentCase, caseFile, err)
result = mtex_base_result(campaign, currentCase, caseFile, "failed");
result.results = struct( ...
    "error_identifier", string(err.identifier), ...
    "error_message", string(err.message));
result.notes = ["MTEX parity case failed during generation."];
end

function result = mtex_base_result(campaign, currentCase, caseFile, status)
result = struct();
result.schema_id = "pytex.parity_result";
result.schema_version = "0.1.0";
result.campaign_id = string(campaign.campaign_id);
result.case_id = string(currentCase.case_id);
result.case_status = string(status);
result.producer = struct( ...
    "system", "mtex", ...
    "system_version", mtex_detect_version(), ...
    "runtime", version, ...
    "script", "scripts/mtex_generators/run_mtex_parity_campaign.m", ...
    "created_utc", string(datetime("now", "TimeZone", "UTC", "Format", "yyyy-MM-dd'T'HH:mm:ss'Z'")));
result.conventions = campaign.conventions;
result.phase = currentCase.phase;
result.tolerances = currentCase.tolerances;
result.results = struct();
result.artifacts = [];
result.provenance = struct( ...
    "input_sha256", mtex_sha256(caseFile), ...
    "case_file", string(caseFile), ...
    "target_baseline", campaign.target_baseline);
result.notes = string.empty(0, 1);
end

function value = get_optional_string(payload, fieldName)
if isfield(payload, fieldName)
    value = string(payload.(fieldName));
else
    value = "";
end
end

