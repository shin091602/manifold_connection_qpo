function [ri,F] = interpolate_torus(fin_qpos,tht0,tht1,tht0_n,tht1_n,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolation of QPOs
%% By:Soi Yamaguchi

%%% input
%fin_qpos :invariant curve solution
%tht0, tht1 :meshgrid
%tht0_n,tht1_n :interpolated grid
%p :parameter

%%% output
%ri :interpolated x,y,z
%F :interpolation function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% torus function matrix
Uinter = fin_qpos(1:p("d")*p("N")*p("M"));
Uinter = reshape(Uinter,p("d")*p("N"),p("M"));
Uinter_full = cell(p("d"),1);
F = cell(p("d"),1);
for dim=1:p("d")
    % dimension index
    i=dim:p("d"):p("N")*p("d");
    % meshgrid
    Uinter_full{dim,1} = repmat(Uinter(i,:)',3,3);
    vi = Uinter_full{dim,1};
    % interpolation
    Fi = griddedInterpolant({tht0,tht1},vi,'spline');
    % store
    F{dim,1} = Fi;
end

ri = cell(p("d"),1);
for dim=1:p("d")
    Fi = F{dim,1};
    ri{dim,1} = Fi({tht0_n,tht1_n});
end

end