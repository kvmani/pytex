function versionText = mtex_detect_version()
%MTEX_DETECT_VERSION Best-effort MTEX version detection across releases.
try
    if exist("mtexVersion", "file")
        versionText = string(mtexVersion);
        return
    end
catch
end
try
    if exist("mtex_version", "file")
        versionText = string(mtex_version);
        return
    end
catch
end
try
    versionText = string(getMTEXpref("version"));
catch
    versionText = "unknown";
end
end

