function locked = adaptive_blend_lock_baseline()
%ADAPTIVE_BLEND_LOCK_BASELINE Load the gateway without creating a session.
% Normally the baseline is unlocked. A local compiler-runtime workaround
% may pin the MEX independently of its balanced per-session locks. Tests
% require returning to this baseline and separately reject closed handles.
% Closing an unused handle is a no-op and loads the gateway without a
% live session. Native stale-handle rejection is checked in the solver test.
surfsmooth3d_adaptive_blend_routs('close',flintmax-1);
locked = mislocked('surfsmooth3d_adaptive_blend_routs');
end
