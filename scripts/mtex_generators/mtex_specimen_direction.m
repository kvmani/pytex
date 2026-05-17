function direction = mtex_specimen_direction(name)
%MTEX_SPECIMEN_DIRECTION Convert a shared specimen direction label to MTEX vector3d.
switch upper(string(name))
    case {"RD", "X"}
        direction = xvector;
    case {"TD", "Y"}
        direction = yvector;
    case {"ND", "Z"}
        direction = zvector;
    otherwise
        error("PyTexParity:UnknownSpecimenDirection", "Unknown specimen direction: %s", name);
end
end

