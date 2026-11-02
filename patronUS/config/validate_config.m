function validate_config(versions, reg)
%VALIDATE_CONFIG Fail early with a readable message on unknown selections.

    required = {'constraints', 'objective'};
    optional = {'dynamics', 'transcription', 'initial_guess'}; % Trajectory segments only
    for k = 1:numel(required)
        if ~isfield(versions, required{k})
            error('validate_config:missing', 'versions.%s is not set.', required{k});
        end
    end

    present = [required, optional(isfield(versions, optional))];
    for k = 1:numel(present)
        f = present{k};
        selected = cellstr(versions.(f)); % A string array lists several (e.g. transcription stages)
        for j = 1:numel(selected)
            key = selected{j};
            if ~isKey(reg.(f), key)
                error('validate_config:unknown', ...
                    'Unknown %s "%s". Available: %s', ...
                    f, key, strjoin(keys(reg.(f)), ', '));
            end
        end
    end

end