function results = cruise_postp_2Wings_dbg(sol_obj,X_struct,dbg,params,fval)

    varNames = fieldnames(X_struct);
    for i = 1:numel(varNames)
        varName = varNames{i};
        results.opt.(varName) = sol_obj.value(X_struct.(varName));
    end

    %% Extract results out of the solution

    results.alpha = sol_obj.value(dbg.alpha);
    results.xi_dw = sol_obj.value(dbg.xi_dw);

    results.J1    = sol_obj.value(dbg.J1);
    results.J2    = sol_obj.value(dbg.J2);
    results.phi1  = sol_obj.value(dbg.phi1);
    results.phi2  = sol_obj.value(dbg.phi2);

    results.CT1   = sol_obj.value(dbg.CT1);
    results.CT2   = sol_obj.value(dbg.CT2);
    results.CP1   = sol_obj.value(dbg.CP1);
    results.CP2   = sol_obj.value(dbg.CP2);
    results.CH1   = sol_obj.value(dbg.CH1);
    results.CH2   = sol_obj.value(dbg.CH2);

    results.CL_wing = sol_obj.value(dbg.CL_wing);
    results.CD_wing = sol_obj.value(dbg.CD_wing);
    results.CM_wing = sol_obj.value(dbg.CM_wing);

    results.CL_tail = sol_obj.value(dbg.CL_tail);
    results.CD_tail = sol_obj.value(dbg.CD_tail);
    results.CM_tail = sol_obj.value(dbg.CM_tail);

    results.CL_fuselage = sol_obj.value(dbg.CL_fuselage);
    results.CD_fuselage = sol_obj.value(dbg.CD_fuselage);
    results.CM_fuselage = sol_obj.value(dbg.CM_fuselage);

    results.L_wing     = sol_obj.value(dbg.L_wing);
    results.L_tail     = sol_obj.value(dbg.L_tail);
    results.L_fuselage = sol_obj.value(dbg.L_fuselage);

    results.D_wing     = sol_obj.value(dbg.D_wing);
    results.D_tail     = sol_obj.value(dbg.D_tail);
    results.D_fuselage = sol_obj.value(dbg.D_fuselage);

    results.MA_wing         = sol_obj.value(dbg.MA_wing);
    results.MA_fuselage     = sol_obj.value(dbg.MA_fuselage);
    results.MA_tail         = sol_obj.value(dbg.MA_tail);

    results.T1 = sol_obj.value(dbg.T1);
    results.T2 = sol_obj.value(dbg.T2);
    results.P1 = sol_obj.value(dbg.P1);
    results.P2 = sol_obj.value(dbg.P2);
    results.H1 = sol_obj.value(dbg.H1);
    results.H2 = sol_obj.value(dbg.H2);

    results.L_Tot = sol_obj.value(dbg.L_Tot);
    results.D_Tot = sol_obj.value(dbg.D_Tot);
    results.W = params.g*params.mass;

    results.x_range = 1/fval *params.prop.diameter* params.prop.max_rps /(params.prop.P_max_eng) * params.prop.eff * params.E_batt;
    results.propeff1  = results.J1*results.CT1/results.CP1 * 100;
    results.propeff2  = results.J2*results.CT2/results.CP2 * 100;

    results.alpha_nac1 = sol_obj.value(dbg.alpha_nac1);
    results.alpha_nac2 = sol_obj.value(dbg.alpha_nac2);
    results.lambda_i1  = sol_obj.value(dbg.lambda_i1);
    results.lambda_i2  = sol_obj.value(dbg.lambda_i2);
    results.CD_nac1    = sol_obj.value(dbg.CD_nac1);
    results.CD_nac2    = sol_obj.value(dbg.CD_nac2);
    results.D_nac1     = sol_obj.value(dbg.D_nac1);
    results.D_nac2     = sol_obj.value(dbg.D_nac2);
    results.vi1        = sol_obj.value(dbg.vi1);
    results.vi2        = sol_obj.value(dbg.vi2);
    results.xi_nac1 = results.phi1-results.alpha_nac1;
    results.xi_nac2 = results.phi2-results.alpha_nac2;
    results.V_nac1 = sqrt(sol_obj.value(dbg.V1_nac2));
    results.V_nac2 = sqrt(sol_obj.value(dbg.V2_nac2));

    results.AeroEff_wing     = results.CL_wing/results.CD_wing;
    results.AeroEff_fuselage = results.CL_fuselage/results.CD_fuselage;
    results.AeroEff_tail     = results.CL_tail/results.CD_tail;

    results.AeroEff_Tot      = results.L_Tot/results.D_Tot;

    results.L_tailwing_wind  = sol_obj.value(dbg.L_tailwing_wind);
    results.D_tailwing_wind  = sol_obj.value(dbg.D_tailwing_wind);

    results.MA               = sol_obj.value(dbg.MA);

    results.M_Fus            = sol_obj.value(dbg.M_Fus);
    results.M_Wing           = sol_obj.value(dbg.M_Wing);
    results.M_Tail           = sol_obj.value(dbg.M_Tail);

    results.M_Eng1           = sol_obj.value(dbg.M_Eng1);
    results.M_Eng2           = sol_obj.value(dbg.M_Eng2);

    results.M_Nac1           = sol_obj.value(dbg.M_Nac1);
    results.M_Nac2           = sol_obj.value(dbg.M_Nac2);

    %% Print results for clarity

    fprintf('\n--- SOLUTION FOUND ---\n');
    fprintf('V =  %10.2f m/s | epsilon1 = %.3f deg | epsilon1 = %.3f deg\n', ...
        results.opt.V, rad2deg(results.opt.epsilon1), rad2deg(results.opt.epsilon2));
    fprintf('alpha = %.3f deg | phi1 = %.3f deg | phi2 = %.3f deg\n', ...
        rad2deg(results.alpha), rad2deg(results.phi1), rad2deg(results.phi2));
    fprintf('n1 = %.3f rps | n2 = %.3f rps\n', ...
        results.opt.n1, results.opt.n2);

    fprintf('\n--- FORCE BREAKDOWN [N] ---\n');
    fprintf('T1 = %10.2f   T2 = %10.2f   H1 = %10.2f  H2 = %10.2f\n', results.T1, results.T2, ...
        results.H1, results.H2);
    fprintf('D_wing = %8.2f  D_fus = %8.2f  D_tail = %8.2f  D_nac1 = %8.2f  D_nac2 = %8.2f  D_Tot = %8.2f\n', ...
        results.D_wing, results.D_fuselage, results.D_tail, ...
        results.D_nac1, results.D_nac2, results.D_Tot);
    fprintf('L_wing = %8.2f  L_fus = %8.2f  L_tail = %8.2f  L_Tot = %8.2f  Weight = %8.2f\n', ...
        results.L_wing, results.L_fuselage, results.L_tail, ...
        results.L_Tot, results.W);

    fprintf('\n--- MOMENT BREAKDOWN [N*m] ---\n');
    fprintf('MA_wing = %8.2f  MA_fus = %8.2f  MA_tail = %8.2f\n', results.MA_wing, ...
        results.MA_fuselage, results.MA_tail);
    fprintf('M_Wing = %8.2f  M_Fus = %8.2f  M_Tail = %8.2f\n', results.M_Wing, results.M_Fus, results.M_Tail);
    fprintf('M_Eng1 = %8.2f  M_Eng2 = %8.2f  M_Nac1 = %8.2f  M_Nac2 = %8.2f\n', ...
        results.M_Eng1, results.M_Eng2, results.M_Nac1, results.M_Nac2);

    fprintf('\n--- CEQ RESIDUALS ---\n');
    fprintf('Long force (should be 0): %10.3f N\n', sol_obj.value(dbg.ceq(1)));
    fprintf('Trans force (should be 0): %10.3f N\n', sol_obj.value(dbg.ceq(2)));
    fprintf('Moment (should be 0):     %10.3f N*m\n', sol_obj.value(dbg.ceq(3)));

end