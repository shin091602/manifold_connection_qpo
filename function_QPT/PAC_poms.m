function Zpo = PAC_poms(z0,s0,vp0,phi0,pacp,p)
%%% input
%z0 :initial guess
%s0 :initial step length size
%phi0 :initial tangent vector
%vp0 :varying parameter dictionary
%F :error vector function
%DF :error vector Jacobian function
%vpf :varying parameter finalization function
%pacp :PAC parameter dictionary
%p :global dictionary

%%% output
%Zpo :updated periodic solution leveraging by PAC

% compute error and error Jacobian
[Fv,vp0] = F_poms(z0,vp0,p);
[DFv,vp0] = DF_poms(z0,vp0,p);

% pac solution package
z = cell(5,2);
z{1,1} = "z";z{1,2} = z0;% solution vector
z{2,1} = "vp0";z{2,2} = vp0;% varying parameter
z{3,1} = "z0";z{3,2} = z0;% previous vector
z{4,1} = "s";z{4,2} = s0;% tangent vector
z{5,1} = "phi0";z{5,2} = phi0;% tangent basis

% save Xd
cd mat/
save("z.mat","z");
cd ..

% correct the initial guess
z = PAC_correction_po(z,pacp,p);

% coflg
if pacp("coflg")
  Zpo = z{1,2};
  return
end

end

