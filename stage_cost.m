
function [lt, Dlx, Dlu, Dlxx, Dluu, Dlux] = stage_cost(xx_t, uu_t, xx_ref, uu_ref, params)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Argument
    %       - state xx at time t
    %       - input uu at time t
    %       - state reference signal xx_ref at t
    %       - input reference signal uu_ref at t
    %   Return
    %       - current cost 
    %       - gradient of cost wrt x, at xx,uu
    %       - gradient of cost wrt u, at xx,uu
    %       - hessian of cost wrt xx, uu, ux
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    QQ = params.cost.QQ;
    RR = params.cost.RR;
    
    [state_dim, ~] = size(xx_ref);
    [input_dim, ~] = size(uu_ref);
    
    % The stage cost is simply defined as a weighted square norm of the
    % errors on the state and input errors 
    lt = (xx_t-xx_ref)'*QQ*(xx_t-xx_ref) + (uu_t-uu_ref)'*RR*(uu_t-uu_ref);

    Dlx = 2*QQ*xx_t - 2*QQ*xx_ref; % Gradient of the cost dL_x = 2Qx
    Dlu = 2*RR*uu_t - 2*RR*uu_ref; % Gradient of the cost dL_u = 2Ru
 
    Dlxx = 2*QQ;
    Dluu = 2*RR;
    Dlux = zeros(input_dim,state_dim);
    
end