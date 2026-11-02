function set_ocp_initial_guess(opti, vars, guess)
%SET_OCP_INITIAL_GUESS Warm start from a previous solution.
%   guess is a results struct from unpack_ocp_solution (e.g. a trapezoidal
%   run). It is interpolated in normalised time tau = t/tf, so the previous
%   run may use a different N or transcription method. Overrides the guess
%   set by define_ocp_variables; names missing in guess keep that guess.

tau_g = guess.t/guess.tf;

for i = 1:vars.nx
    name = vars.state_names{i};
    if isfield(guess, name)
        opti.set_initial(vars.X(i,:), interp1(tau_g, guess.(name), vars.tau));
    end
end

for i = 1:vars.nu
    name = vars.control_names{i};
    if isfield(guess, name)
        opti.set_initial(vars.U(i,:), interp1(tau_g, guess.(name), vars.tau));
        if size(vars.Um, 2) > 0
            opti.set_initial(vars.Um(i,:), interp1(tau_g, guess.(name), vars.tau_mid));
        end
    end
end

opti.set_initial(vars.tf, guess.tf);

end
