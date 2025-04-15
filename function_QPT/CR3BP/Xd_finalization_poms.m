function Xd = Xd_finalization_poms(z_mul,Xd,p)
%%% input
% z_mul(1:d*M) :X of patch points
% z_mul(d*M+1) :the period of the periodic orbit
% Xd: data package:xT,phiT,f
% p:variables in dictionary

%%% output
% Xd: updated data package:xT,phiT,f,x0,f0

%% DICTIONARY
d = p("d");
mu = p("mu");

%% UPDATE PACKAGE
% update x0
x0 = z_mul(1:d);
Xd{1,2} = x0;
% update f0
f0 = fun_cr3bp([],x0,mu);
Xd{2,2} = f0;

end