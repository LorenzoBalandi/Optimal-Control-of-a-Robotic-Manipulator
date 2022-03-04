

function [lT, Dlx, Dlxx] = term_cost(xx_t, xx_ref, params)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Argument
    %       - state xx at time T
    %       - reference signal xx_ref at T
    %   Return
    %       - current cost 
    %       - gradient of cost wrt x, at xx_T
    %       - Hessian of the term cost
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     gg = params.dyn.gg;
%     ll = params.dyn.ll; 
%     kk = params.dyn.kk;
%     mm = params.dyn.mm;
%     
%     dt = params.dyn.dt;
%     
%     [state_dim, ~] = size(xx_ref);
    
%      QQ = params.cost.QQ;
%      RR = params.cost.RR;
%      pp = zeros(4,1);
%      uu_ref = [24.956;4.1459];
%      [~, fx_cost, fu_cost, pfxx_dummy, pfuu_dummy, pfux_dummy] = dynamics(xx_ref(:,1), uu_ref(:,1), pp);
%      AA = fx_cost'; 
%      BB = fu_cost';
%      QQf = idare(AA,BB,QQ,RR);

    QQf = params.cost.QQf;

    lT = 0.5*(xx_t-xx_ref)'*QQf*(xx_t-xx_ref);

    Dlx = QQf*xx_t - QQf*xx_ref; % Gradient of the cost: dL_x = 2Qx
    
    Dlxx = QQf;

end