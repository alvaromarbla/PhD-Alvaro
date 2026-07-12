function [models] = models_config(versions)

if versions.AC_version == "1WING"

models.CL_lookup = cl_model(versions.CL_version);
models.CD_lookup = cd_model(versions.CD_version);
models.CT_lookup = ct_model(versions.CT_version);
models.CP_lookup = cp_model(versions.CP_version);
models.CL_fuselage_lookup = cl_fuselage(versions.CL_version);
models.CD_fuselage_lookup = cd_fuselage(versions.CD_version);

elseif AC_version == "2WING"

models.CL_lookup = cl_wing(versions.CL_version);
models.CD_lookup = cd_wing(versions.CD_version);
models.CT_lookup = ct_model(versions.CT_version);
models.CP_lookup = cp_model(versions.CP_version);

models.CL_fuselage_lookup = cl_fuselage(CL_version);
models.CD_fuselage_lookup = cd_fuselage(versions.CD_version);
else
    warning('AC version unknown')

end