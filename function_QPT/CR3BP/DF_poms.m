function [DF_poms,Xd] = DF_poms(z_mul,Xd,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jacobian for periodic orbit using mutiple shooting
%% By:Soi Yamaguchi
%%% input
% z_mul(1:d*M) :X of patch points
% z_mul(d*M+1) :the period of the periodic orbit
% Xd: data package:x0,f0,xT,phiT,f
% p:parameter dictionary

%%% output
% DF_poms:Jacobian of Error vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY
d = p("d");
M = p("M");
mu = p("mu");

%% PACKAGE
phiT = Xd{4,2};
f = Xd{5,2};

%% ERROR PARTIAL
DF = zeros(d*M,d*M+1);
for i=2:M
  j=i-1;%j:1~M-1
  % partial of integrated patch points
  DF(d*j-(d-1):d*j,d*j-(d-1):d*j) = phiT(:,:,j);
  % partial of patch point
  DF(d*j-(d-1):d*j,d*i-(d-1):d*i) = -eye(d,d);
  % partial of period
  DF(d*j-(d-1):d*j,d*M+1) = f(:,j)/M;% divided by ...patch to patch time:T/M
end

% in case of last patch point to first patch point
% partial of integrated patch points
DF(d*M-(d-1):d*M,d*M-(d-1):d*M) = phiT(:,:,end);
% partial of patch point
DF(d*M-(d-1):d*M,1:d) = -eye(d,d);
% partial of period
DF(d*M-(d-1):d*M,d*M+1) = f(:,end)/M;

% phase condition partials
if p("phsflg")
  Dphs = zeros(1,d*M+1);
  Dphs(1,1:d) = Xd{2,2}';
end

% parametrization contraints partials
% Jacobi Integral partials
DCcon_po = [Jacobi_const_diff(z_mul,mu)',zeros(1,(M-1)*d+1)];

% append
DF_poms = [DF;Dphs;DCcon_po];

end
