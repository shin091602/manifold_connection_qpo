function z = PAC_correction_qpo(z,zpack,zpo,pacqp,p,FR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correction via a Newton's method to calculate QPOs
%% By:Soi Yamaguchi

%%% input
%z :initial guess
%zpack :varying parameter for finalization
%s0 :initial step length size
%phi0 :initial tangent vector
%vp0(Xd) :varying parameter
%F :error vector function
%DF :error vector Jacobian function
%vpf :varying parameter finalization function
%pacqp :PAC parameter
%FR :fourier matrix data package

%%% output
%z :updated periodic solution leveraging by correction calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Error vector
[Fval,Xd] = F_qpoms(z{1,2},zpo,z{2,2},p,FR);
% update fT,Uf --final state
z{2,2} = Xd;

% add pseudo-arclength continuation
% use current and previous solution
% (Ref:(p172 s1_z))
s = (1/pacqp("N")).*dot((z{1,2}(1:pacqp("fpidx")))-z{3,2}(1:pacqp("fpidx")),z{5,2}(1:pacqp("fpidx")))+dot(z{1,2}(pacqp("fpidx")+1:end)-z{3,2}(pacqp("fpidx")+1:end),z{5,2}(pacqp("fpidx")+1:end))-z{4,2};
% renew error vector
%(Ref:(23))
Fval_n = [Fval;s];

% correct guess
for iteration=1:pacqp("itmax")
    % compute Jacobian
    [DFval,Xd] = DF_qpoms(z{1,2},z{2,2},p,FR);
    z{2,2} = Xd;

    % add pseudo-arclength continuation Jacobian
    Ds = [(1/pacqp("N")).*(z{5,2}(1:pacqp("fpidx")))',z{5,2}(pacqp("fpidx")+1:end)'];
    % renew error vector
    DFval_n = [DFval;Ds];

    %correction
    dz = sparse(-DFval_n)\Fval_n;

    %apply correction
    z{1,2} = z{1,2}+dz;

    %updated error
    [Fval,Xd] = F_qpoms(z{1,2},zpo,z{2,2},p,FR);
    z{2,2} = Xd;

    % add pseudo-arclength continuation
    % use current and previous solution
    s = (1/pacqp("N")).*dot((z{1,2}(1:pacqp("fpidx")))-z{3,2}(1:pacqp("fpidx")),z{5,2}(1:pacqp("fpidx")))+dot(z{1,2}(pacqp("fpidx")+1:end)-z{3,2}(pacqp("fpidx")+1:end),z{5,2}(pacqp("fpidx")+1:end))-z{4,2};
    % renew error vector
    Fval_n = [Fval;s];
    err = norm(Fval_n);
    disp(err)
    if err<pacqp("tol")
        break;
    end
end

if (err>pacqp("tol"))&&(iteration==pacqp("itmax"))
    disp("Solution not found");
    z{1,2} = [];
    return
else
    disp('PAC solution converged')
    disp(strcat('Number of iterations:',num2str(iteration)))
    disp(strcat('Final Error:',num2str(err)))
    disp(strcat('Step Size:',num2str(z{4,2})))

    [DFval,Xd] = DF_qpoms(z{1,2},z{2,2},p,FR);
    z{2,2} = Xd;
    Ds = [(1/pacqp("N")).*(z{5,2}(1:pacqp("fpidx")))',z{5,2}(pacqp("fpidx")+1:end)'];
    % renew error vector
    DFval_n = [DFval;Ds];
    % nullspace
    % update tangent basis
    phi = sparse(DFval_n)\[zeros(length(DFval_n)-1,1);1];
    z{5,2} = phi/sqrt((dot(phi(1:pacqp("fpidx")),phi(1:pacqp("fpidx")))+dot(phi(pacqp("fpidx")+1:end),phi(pacqp("fpidx")+1:end))));

    % update step length
    Eps = abs(pacqp("optit")/iteration);
    if (Eps>2)||(iteration==0)
        Eps = 2;
    else
        if Eps<0.5
            Eps = 0.5;
        end
    end
    if z{4,2}>0
        z{4,2} = min(pacqp("smax"),z{4,2}*Eps);
    else
        z{4,2} = max(pacqp("smax"),z{4,2}*Eps);
    end
end
z{2,2} = Xd_finalization_qpoms(z{1,2},zpack,p,FR);

end
