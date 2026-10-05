function label = sigma_mode_label(flag)
%SIGMA_MODE_LABEL Readable names for the existing sigma option path.
switch flag
    case 0
        label = 'Constant (non-adaptive)';
    case 1
        label = 'Adaptive: longest triangle side';
    case 2
        label = 'Adaptive: recursive (longest-side seeds)';
    case 3
        label = 'Adaptive: shortest triangle side';
    otherwise
        error('edgepreserve:sigmaMode','Unknown sigma mode: %g',flag);
end
end
