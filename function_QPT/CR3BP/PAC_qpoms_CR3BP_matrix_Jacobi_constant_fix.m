function sol_qpos = PAC_qpoms_CR3BP_matrix_Jacobi_constant_fix(z0,zpo,s0,vp0,phi0,p,pacqp, C_periodic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pseudo-arclength continuation and correction in computing quasi-periodic orbits

%%% input
%z0 :initial guess
%zpo :periodic solution
%s0 :initial step length size
%phi0 :tangent vector
%pacqp :paramters

%%% output
%sol_qpos :quasi-peiodic solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute error vector and jacobian
[~,vp0] = F_qpoms_CR3BP_matrix_Jacobi_constant_fix(z0,zpo,vp0,p,C_periodic);
[~,vp0] = DF_qpoms_CR3BP_matrix_Jacobi_constant_fix(z0,vp0,p,zpo);

% PAC data package
ps = cell(5,2);
ps{1,1} = "z";ps{1,2} = z0;% solution vector
ps{2,1} = "vp0";ps{2,2} = vp0;% varying parameter
ps{3,1} = "z0";ps{3,2} = z0;% previous vector
ps{4,1} = "s";ps{4,2} = s0;% tangent vector
ps{5,1} = "phi0";ps{5,2} = phi0;% tangent basis

%% CORRECTION
% correct the initial guess
ps = PAC_correction_qpo_CR3BP_matrix_Jacobi_constant_fix(ps,vp0,zpo,pacqp,p,C_periodic);
disp("corrected")

%% CONTINUATION
[DFv,~] = DF_qpoms_CR3BP_matrix_Jacobi_constant_fix(ps{1,2},ps{2,2},p,zpo);
% nullspace
SVD_null = svd(DFv);
null_acuracy = min(SVD_null)*1e+1;
ps{5,2} = null(DFv,null_acuracy);
null_disp = strcat("Null size:",num2str((size(ps{5,2}))));
disp(null_disp)
while (size(ps{5,2})) > 1
    null_acuracy = null_acuracy*0.5;
    ps{5,2} = null(DFv,null_acuracy);
end
null_disp = strcat("Null size:",num2str((size(ps{5,2}))));
disp(null_disp)

% arrange the angle of vector in solution space
if length(phi0)==length(ps{5,2})
    if acos(dot(ps{5,2},phi0)/(norm(ps{5,2})*norm(phi0)))>pi/2
        ps{5,2} = -ps{5,2};
    end
end

% store solution
sol_qpos = cell(pacqp("n"),1);

% initial solution
solc_qpos = cell(5,2);
solc_qpos{1,1} = "z";solc_qpos{1,2} = ps{1,2};% solution vector
solc_qpos{2,1} = "vp0";solc_qpos{2,2} = ps{2,2};% varying parameter
solc_qpos{3,1} = "z0";solc_qpos{3,2} = ps{3,2};% previous vector
solc_qpos{4,1} = "s";solc_qpos{4,2} = ps{4,2};% tangent vector
solc_qpos{5,1} = "phi0";solc_qpos{5,2} = ps{5,2};% tangent basis
sol_qpos{1,1} = solc_qpos;

% predict next solution by PAC
for i=2:pacqp("n")
    % previous solution
    ps{3,2} = ps{1,2};
    % continuation
    ps{1,2} = ps{3,2}+ps{5,2}*ps{4,2};
    % correct solution
    ps = PAC_correction_qpo_CR3BP_matrix_Jacobi_constant_fix(ps,ps{2,2},zpo,pacqp,p,C_periodic);

    % break
    if isempty(ps{1,2})==1
        sol_qpos = sol_qpos{1:i-1,1};
        break
    end

    % store current solution
    solc_qpos = cell(5,2);
    solc_qpos{1,1} = "z";solc_qpos{1,2} = ps{1,2};% solution vector
    solc_qpos{2,1} = "vp0";solc_qpos{2,2} = ps{2,2};% varying parameter
    solc_qpos{3,1} = "z0";solc_qpos{3,2} = ps{3,2};% previous vector
    solc_qpos{4,1} = "s";solc_qpos{4,2} = ps{4,2};% tangent vector
    solc_qpos{5,1} = "phi0";solc_qpos{5,2} = ps{5,2};% tangent basis
    sol_qpos{i,1} = solc_qpos;

end

end
