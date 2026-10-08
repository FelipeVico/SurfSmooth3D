function enabled = validate_newton_recovery(value)
%VALIDATE_NEWTON_RECOVERY Validate the public guarded-solver switch.
if (~isnumeric(value) && ~islogical(value)) || ~isscalar(value) || ...
        ~isreal(value) || ~ismember(value,[0 1])
    error('MULTISCALE_MESHER:InvalidNewtonRecovery', ...
        'newton_recovery must be a scalar logical or numeric zero/one.');
end
enabled = logical(value);
end
