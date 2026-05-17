function payload = mtex_read_json(path)
%MTEX_READ_JSON Read a UTF-8 JSON file into a MATLAB struct.
text = fileread(path);
payload = jsondecode(text);
end

