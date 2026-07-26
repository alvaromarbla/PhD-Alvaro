function validate_config(versions, reg)
%VALIDATE_CONFIG Fail early with a readable message on unknown selections.

    required = {'constraints', 'objective'};
    for k = 1:numel(required)
        f = required{k};
        if ~isfield(versions, f)
            error('validate_config:missing', 'versions.%s is not set.', f);
        end
        key = char(versions.(f));
        if ~isKey(reg.(f), key)
            error('validate_config:unknown', ...
                'Unknown %s "%s". Available: %s', ...
                f, key, strjoin(keys(reg.(f)), ', '));
        end
    end

end