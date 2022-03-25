%% Safety Controller matlab version for evaluation pyton code

a = 2;
n = 200;
time_step = 0.0001; % 0.1ms
% position, velocity and acceleration column vectors
x_t_wo = 0;x_t_wo1 = 0;x_t_wo2 = 0;
v_t_wo = 0;v_t_wo1 = 0;v_t_wo2 = 0;
a_t_wo = 0;a_t_wo1 = 0;a_t_wo2 = 0;


% Safety state
x_t_s = 0;x_t_s_1 = 0;x_t_s_2 = 0;
v_t_s = 0;v_t_s_1 = 0;v_t_s_2 = 0;
a_t_s = 0;a_t_s_1 = 0;a_t_s_2 = 0;

x_eq = 0;
%
b_h = 0.0; % human damping variable
gamma_safe = 10;
% Constraint parameters
max_joint_acceleration = 100;
max_joint_limit = 1.57;
max_b_r = 100;
%Impedance variable
I = 0.01;
K = 0;
B_OLD = 0.1;
for i=1:200
    %Input torque
    torque_applied = a*sin(i*(2*pi/n));
    
    save_torque_applied(i) = torque_applied;
    %Without safety just x axis
    save_xt(i) = x_t_wo;
    save_vt(i) = v_t_wo;
    save_at(i) = a_t_wo;
    %After safety just x axis
    save_xt_s(i) = x_t_s;
    save_vt_s(i) = v_t_s;
    save_at_s(i) = a_t_s;
    
    %% Safety Controller
    % Before Safetycontroller
    [x_t, v_t, a_t] = Doadmittance(time_step, I, B_OLD, K, x_eq, x_t_s_1, x_t_s_2,v_t_s_1,torque_applied);
    % Safety Controller  
    [F,Lambda,B_mat,E_mat] = calControlAffineSysMatrix(x_t,v_t,I,b_h,K,x_eq, torque_applied);
    [h,dh_dx,L_fh,L_Fhh,L_gh,L_haf,L_ghu1] = MakeCbf(x_t,v_t,F,Lambda,B_mat,E_mat,max_joint_limit,max_joint_acceleration);
    [sol] = Qpsolver(x_t,v_t,a_t,B_OLD,F,Lambda,B_mat,E_mat,h,dh_dx,L_fh,L_Fhh,L_gh,L_haf,L_ghu1,max_joint_acceleration,max_joint_limit,max_b_r,b_h,gamma_safe);
    B_SAFE = sol(1); %% sol(1) our new control input
    
    % After Safetycontroller
    [x_t_s, v_t_s, a_t_s] = SafeDoadmittance(time_step, I, B_SAFE, K, x_eq, x_t_s_1, x_t_s_2,v_t_s_1,torque_applied,sol);
    
    %% Without Safety Controller
    [x_t_wo, v_t_wo, a_t_wo] = Doadmittance(time_step, I, B_OLD, K, x_eq, x_t_wo1, x_t_wo2,v_t_wo1,torque_applied);
    
    %% SAVE the WITHOUT SAFETY sate 
    % save t-1 value to t-2 vlaue
    x_t_wo2 = x_t_wo1;
    v_t_wo2 = v_t_wo1;
    a_t_wo2 = a_t_wo1;
    % save t value to t-1 value
    x_t_wo1 = x_t_wo;
    v_t_wo1 = v_t_wo;
    a_t_wo1 = a_t_wo;
    
    %% SAVE the WITH SAFETY states
    x_t_s_2 = x_t_s_1;
    v_t_s_2 = v_t_s_1;
    a_t_s_2 = a_t_s_1;
    % save t value to t-1 value
    x_t_s_1 = x_t_s;
    v_t_s_1 = v_t_s;
    a_t_s_1 = a_t_s;
    
    save_b_old(i) = B_OLD;
    save_b_new(i) = B_SAFE;
    
end


% subplot(1,1,4)
% subplot(1,2,4)
% subplot(2,1,4)
% subplot(2,2,4)

function [x_t, v_t, a_t] = Doadmittance(time_step, I, B, K, x_eq, x_t_1, x_t_2,v_t_1,force_applied)
        force_applied = -force_applied;
        t = time_step;
        A_dyn =  (I / (t * t) + B / (t) + K);
        B_temp = I*((2.0 * x_t_1 - x_t_2) / (t * t));
        C_temp = B*(x_t_1 / (t));
        D_temp = K*x_eq;
        
        x_t = A_dyn\(-force_applied + B_temp + C_temp + D_temp);
        
        v_t = (x_t-x_t_1)/(t);
        
        a_t = (v_t - v_t_1)/(t);       
end

function [x_t_s, v_t_s, a_t_s] = SafeDoadmittance(time_step, I, B, K, x_eq, x_t_1, x_t_2,v_t_1,force_applied,sol)
        force_applied = -force_applied;
        t = time_step;
        u1 = sol(3);
        A_dyn =  (I / (t * t) + B / (t) + K);
        B_temp = I*((2.0 * x_t_1 - x_t_2) / (t * t));
        C_temp = B*(x_t_1 / (t));
        D_temp = K*x_eq;
        
        x_t_s = A_dyn\(-force_applied + u1 + B_temp + C_temp + D_temp);
        
        v_t_s = (x_t_s - x_t_1)/(t);
        
        a_t_s = (v_t_s - v_t_1)/(t);
           
end

function [F,Lambda,B_mat,E_mat] = calControlAffineSysMatrix(x,v,I,b_h,K,x_eq, force_applied) 
        
    F = [0 1;
        -I\K 0];

    Lambda = [0;
        -I\v]*b_h;

    B_mat = [0;
           -I\v];

    E_mat = [0;
        I\(force_applied+K*x_eq)];   % Before this part have - sign so I changed the sing - to + , after changing it works
end


function [h,dh_dx,L_fh,L_Fhh,L_gh,L_haf,L_ghu1] = MakeCbf(x,v,F,Lambda,B_mat,E_mat,max_joint_limit,max_joint_acceleration)
    % this function make reciprocal control barrier function h(x)>=0
    % h=(y_m-sign(x(2))*x(1))-((1/2)*(x(2)^2)/(a_m));


    h = max_joint_limit - sign(v)*x - 1/2*(v^2/max_joint_acceleration);

    

    dh_dx = [-sign(v);
            -v/max_joint_acceleration];
    
    %% CBF Constraints  -dh/dx*g(x)*b_r =< alpha*h + dh/dx*f(x)
    %hdot = dh/dx*xdot = dh/dx(f(x) + g(x)*u) (f(x) : L_fh+L_Fhh+L_haf, g(x) = L_gh)
    %L_fh = dh_dx * (A*x);
    L_fh = dh_dx'*(F*[x;v]);
    
    %L_Fhh  
    L_Fhh = dh_dx'*Lambda;
    
    % L_haf
    L_haf = dh_dx'*E_mat;
    
    % This term g(x) vector lie derivative L_gh = dh_dx*B_mat;
    L_gh = dh_dx'*B_mat;
    
    L_ghu1 = dh_dx'*[0;1];
 

end

function [sol] = Qpsolver(x,v,a,B,F,Lambda,B_mat,E_mat,h,dh_dx,L_fh,L_Fhh,L_gh,L_haf,L_ghu1,max_joint_acceleration,max_joint_limit,max_b_r,b_h,gamma_safe)
    %Just see x axis
    % Decision variable order = control input, slack variable, u1
    HH=zeros(3,3);HH(1,1)=1;HH(2,2)=0;HH(3,3)=1000;
    HH = 2*HH;
    FF=[-B; 10;0];
    %For acceleration constraint
    gg = [0,1;
          0,-1];
    
    
    
    f_x =  F*[x;v] + Lambda + E_mat;  %without g(x)*u term x_dot = f(x) + g(x)*u
    ac_f = gg*f_x;
    ac_g = gg*B_mat;
    ac_g1 = gg*[0;1];
    
    %For CBF constraint
    %h_dot >=0   ->   dh/dx*[f(x) + g(x)*b_r] >= -alpha*h   ->   -dh/dx*g(x)*b_r =< alpha*h + dh/dx*f(x)
    h_temp_x = L_fh + L_Fhh + L_haf + gamma_safe*h;  %%Just see x axis value == alpha*h + dh/dx*f(x)
    
    AA = [0,-1,0;   %% slack varialbe delta>0
          1,-1,0;   %% control input constraint
          -1,-1,0;  %% control input constraint
          -L_gh,0,-L_ghu1; %% CBF constraint his term g(x) vector lie derivative L_gh = dh_dx*B_mat;
          ac_g(1),0,ac_g1(1);   %% accleration constraint
          ac_g(2),0,ac_g1(2)];  %% accleration constraint

    bb = [0;max_b_r;max_b_r;h_temp_x; -ac_f(1) + max_joint_acceleration; -ac_f(2) + max_joint_acceleration];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  HH = 3 by 3, FF = 3 by 1, AA = 6 by 3, bb = 6 by 1   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    options = optimset('Display', 'off');
    [uu,J,exitflag]=quadprog(HH,FF,AA,bb,[],[],[],[],[],options);
    sol = uu
    exitflag
end