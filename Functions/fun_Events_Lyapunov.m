function [value, isterminal, direction] = fun_Events_Lyapunov(t, x, varargin)
% Stop integration when the value of 'value' becomes zero from the direction of 'direction'

    value = x(2);    % The value that we want to be zero
    isterminal = 1;  % Halt integration
    direction = 0;   % The zero can be approached from either direction
    % direction = +1;   % Positive direction only
    % direction = -1;   % Negative direction only
end
