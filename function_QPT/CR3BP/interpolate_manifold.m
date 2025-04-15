function [rim,Fm,tm,Xs_full] = interpolate_manifold(fin_qpos,del_w_us_inter_full,tht0,tht1,tht0_n,tht1_n,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolation of unstable and stable manifolds emataing from QPO
% By:Soi Yamaguchi

%%% input
%fin_qpos :invariant curve solution
%del_w_us_inter_full :unstable direction of manifolds
%tht0, tht1 :meshgrid
%p :parameter

%%% output
%rim :interpolated coordinates x,y,z
%Fm :interpolation function
%tm :transfer time history
%Xs :torus manifolds grid (no interpolate, row data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dictionary
d = p("d");
N = p("N");
M = p("M");
mu = p("mu");

%% OPTIONS ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13);

snap_time = linspace(p("snap_ini_time"),p("snap_fin_time"),p("snap_span"));
% initialization
Xs_full = zeros(d*N*M,length(snap_time));
xini_opt_n = [];
% calculate manifold
    for i=1:M
        U0 = zeros(d*N,1);
        % invariant curve idx
        k=d*N*(i-1)+1:d*N*i;
        X0 = fin_qpos(k);
        del_w_us = del_w_us_inter_full(:,:,i);
        for j=1:N
            % % Create perturbation vector
            % vpert = p("xpert")*norm(X0(d*j-2:d*j))/norm(X0(d*j-5:d*j-3));
            % pert = [ones(3,1).*p("xpert"); ones(3,1).*vpert];
            U0(d*j-5:d*j) = X0(d*j-5:d*j)+p("pert").*del_w_us(:,j);
        end
        % torus manifold grid
        xini_opt_n = [xini_opt_n;U0];
    end
    for sn = 1:length(snap_time)
        for i=1:M
            Xs = zeros(d*N,1);
            % invariant curve idx
            k=d*N*(i-1)+1:d*N*i;
            U0 = xini_opt_n(k);
            for j=1:N
            % time scale
            tspan = linspace(0,snap_time(sn),3);
            [~,X] = ode113(@(t,x) fun_cr3bp(t,x,mu),tspan,U0(d*j-5:d*j),options_ODE);
            Xs(d*j-5:d*j) = X(end,:)';
            end
            % torus manifold grid
            Xs_full(k,sn) = Xs;
        end
    end

tm = [];

% initialization of interpolate function
Fm = cell(d,length(snap_time));
for sn = 1:length(snap_time)
    % torus manifold of each span
    Xs = Xs_full(:,sn);
    Xs = reshape(Xs,d*N,M);
    for dim=1:d
        % dimention index
        i=dim:d:N*d;
        % meshgrid
        vsi = repmat(Xs(i,:)',3,3);
        % interpolation
        Fi = griddedInterpolant({tht0,tht1},vsi,'spline');
        % store
        Fm{dim,sn} = Fi;
    end
end

% interpolation
rim = cell(3,p("snap_span"));
for sn=1:p("snap_span")
    for dim=1:3
        Fi = Fm{dim,sn};
        rim{dim,sn} = Fi({tht0_n,tht1_n});
    end
end

end