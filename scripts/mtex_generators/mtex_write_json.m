function mtex_write_json(path, payload)
%MTEX_WRITE_JSON Write a MATLAB struct as pretty JSON when supported.
try
    text = jsonencode(payload, "PrettyPrint", true);
catch
    text = jsonencode(payload);
end
fid = fopen(path, "w", "n", "UTF-8");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", text);
end

