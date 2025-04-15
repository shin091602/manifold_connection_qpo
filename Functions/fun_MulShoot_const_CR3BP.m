% Constraint
function [c, ceq] = fun_MulShoot_const_CR3BP(X, r0, rf, n, ToF, mu,options)
% X       : design variable vector
% r0      : position of the initial point
% rf      : position of the final point
% n       : the number of the patch points
% ToF     : Time of Flight from the initial to final points
% mu      : mass ratio of the primaries
% options : options for ode

  if(size(r0,2) > 1), r0 = r0'; end
  if(size(rf,2) > 1), rf = rf'; end

  % Equality constraint
  ceq = [sum(X((end-n+2):end)) - ToF];

  % Inequality constraint
  c = [];
end


