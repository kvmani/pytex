function digest = mtex_sha256(path)
%MTEX_SHA256 Return the SHA-256 digest of a file.
engine = java.security.MessageDigest.getInstance("SHA-256");
fid = fopen(path, "r");
cleanup = onCleanup(@() fclose(fid));
while ~feof(fid)
    chunk = fread(fid, 8192, "uint8");
    engine.update(chunk);
end
hash = typecast(engine.digest(), "uint8");
digest = lower(reshape(dec2hex(hash)', 1, []));
end

