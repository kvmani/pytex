function run_mtex_parity_campaign(caseFile, outputRoot)
%RUN_MTEX_PARITY_CAMPAIGN Generate MTEX parity JSON artifacts from a case campaign.
%
%   run_mtex_parity_campaign("fixtures/mtex_parity/campaigns/orientation_core_cases.json", ...
%                            "fixtures/mtex_parity/results/mtex")
%
% This script is intended to run on a separate MATLAB installation with MTEX 6.0.0.
% It reads the same campaign JSON consumed by PyTex, computes MTEX outputs, and writes one
% pytex.parity_result JSON file per active case. Pending cases are emitted as skipped JSON records
% so provenance is explicit.

if nargin ~= 2
    error("PyTexParity:InvalidArguments", "Expected caseFile and outputRoot arguments.");
end

caseFile = string(caseFile);
outputRoot = string(outputRoot);

campaign = mtex_read_json(caseFile);
campaignOutput = fullfile(outputRoot, string(campaign.campaign_id));
if ~exist(campaignOutput, "dir")
    mkdir(campaignOutput);
end

cases = campaign.cases;
caseCount = numel(cases);
for idx = 1:caseCount
    if iscell(cases)
        currentCase = cases{idx};
    else
        currentCase = cases(idx);
    end
    try
        result = mtex_run_single_parity_case(campaign, currentCase, caseFile);
    catch err
        result = mtex_failed_result(campaign, currentCase, caseFile, err);
    end
    resultFile = fullfile(campaignOutput, string(currentCase.case_id) + ".json");
    mtex_write_json(resultFile, result);
end
end
