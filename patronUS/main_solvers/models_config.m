function [models] = models_config(versions)

if versions.AC_version == "1Wing"

models.CL_lookup = cl_model(versions.CL_version);
models.CD_lookup = cd_model(versions.CD_version);
models.CT_lookup = ct_model(versions.CT_version);
models.CP_lookup = cp_model(versions.CP_version);
models.CL_fuselage_lookup = cl_fuselage(versions.CL_version);
models.CD_fuselage_lookup = cd_fuselage(versions.CD_version);

elseif versions.AC_version == "2Wings"

% Aero
models.CL_wing_lookup = cl_wing(versions.CL_version);
models.CD_wing_lookup = cd_wing(versions.CD_version);
models.CM_wing_lookup = cm_wing(versions.CM_version);

models.CL_fuselage_lookup = cl_fuselage(versions.CL_version);
models.CD_fuselage_lookup = cd_fuselage(versions.CD_version);
models.CM_fuselage_lookup = cm_fuselage(versions.CD_version);

models.CL_tail_lookup = cl_tail(versions.CL_version);
models.CD_tail_lookup = cd_tail(versions.CD_version);
models.CM_tail_lookup = cm_tail(versions.CD_version);

% Prop
models.CT_lookup = ct_model(versions.CT_version);
models.CP_lookup = cp_model(versions.CP_version);
models.CH_lookup = ch_model(versions.CP_version);

% Nac

models.CD_nac_lookup = cd_fuselage(versions.CD_version);

else
    warning('AC version unknown')

end