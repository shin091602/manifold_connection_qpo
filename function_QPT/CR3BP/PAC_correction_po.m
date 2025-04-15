function z = PAC_correction_po(z,pacp,p)
%%% input
%z0 :initial guess
%s0 :initial step length size
%phi0 :initial tangent vector
%vp0(Xd) :varying parameter
%F :error vector function
%DF :error vector Jacobian function
%vpf :varying parameter finalization function
%p :PAC parameter

%%% output
%z :updated periodic solution leveraging by correction calculation

%define error vector
[Fval,Xd] = F_poms(z{1,2},z{2,2},p);
z{2,2} = Xd;

% correct guess
for iteration=1:pacp("itmax")
  % compute Jacobian
  [DFval,Xd] = DF_poms(z{1,2},z{2,2},p);
  z{2,2} = Xd;

  %correction
  if pacp("issparse")
    dz = sparse(-DFval)\Fval;
  else
    dz = -DFval\Fval;
  end

  %apply correction
  z{1,2} = z{1,2}+dz;

  %updated error
  [Fval,Xd] = F_poms(z{1,2},z{2,2},p);
  z{2,2} = Xd;

  err = norm(Fval);
  disp(err)
end

if err>pacp("tol")
  disp("Solution not found");
  z{1,2} = [];
end

if err<pacp("tol")
  disp("PAC solution converged.")
  Xd = Xd_finalization_poms(z{1,2},Xd,p);
  z{2,2} = Xd;
  return;
else
  if iteration==pacp("itmax")
    disp("PAC solution Not converged.")
  end
end
end
