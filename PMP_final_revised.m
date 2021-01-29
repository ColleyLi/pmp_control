% PMP Code written by Chonghyuk Song (Date: July 6, 2017)
% Email: chonghyuk.song@kaist.ac.kr

%% 1. Loading All Variables

% Battery
load('Battery.mat');
Q_cap = Bat.Q_cap;
SOC_ind = Bat.SOC_ind;
Res_dis = Bat.Res_dis;
Res_cha = Bat.Res_cha;
Vol_dis = Bat.Vol_dis;          %same as Vol_cha
Cur_lim_dis = Bat.Cur_lim_dis;  %current limit (discharge)
Cur_lim_cha = Bat.Cur_lim_cha;  %current limit (charge)
Pow_max = Bat.Pow_max;
Pow_min = Bat.Pow_min;

% Motor 1
load('Motor1.mat');
sp1 = Mot.sp;
tq1 = Mot.tq;
sp_full1 = Mot.sp_full;
%we need to curb the speed region where the operating limit is drawn
sp_full1(1) = -628.3;
sp_full1(end) = 628.3;

tq_max1 = Mot.tq_max;
tq_min1 = Mot.tq_min;
eff1 = Mot.eff;

% Motor 2
load('Motor2.mat');         % Motor2 is the smaller of the two motors
sp2 = Mot.sp;
tq2 = Mot.tq;
sp_full2 = Mot.sp_full;
tq_max2 = Mot.tq_max;
tq_min2 = Mot.tq_min;
eff2 = Mot.eff;

% Engine
load('Engine.mat');
spE = Eng.sp;
% we would like to extend the feasible speed range of the engine from min
% 105 to minimum 0:
spE = [0:13.0875:13.0875*7, 104.6, spE];

tqE = Eng.tq;
fc = Eng.fc;
%accordingly, we have to add extra information to the fuel consumption map
%as a result of extending the speed range of the engine
fc = [zeros(length(tqE),9), fc];

sp_fullE = Eng.sp_full;
sp_fullE = [0:13.0875:13.0875*7, 104.6, sp_fullE];
tq_maxE = Eng.tq_max;
% the engine torque limit at engine speeds below 105 rad/s is 0
% i.e. the engine is not turned on at below this speed
tq_maxE = [zeros(1,9), tq_maxE];


% Vehicle Specifications
load('Vehicle.mat');
m_veh = VehData.comp.Veh.m_veh;               % (kg)
CD = VehData.comp.Veh.CoeffDrag;              % drag coefficient (dimensionless)
Rtire = VehData.comp.Veh.r_tire;              % (m)
Rr1 = VehData.comp.Veh.rollingRes1;
Rr2 = VehData.comp.Veh.rollingRes2;
Rr3 = VehData.comp.Veh.rollingRes3;
FA  = VehData.comp.Veh.AreaFront;             % (m^2)
SR  = VehData.comp.TM.gearratio;
Gear_eff = VehData.comp.TM.geareff;
r_final = VehData.comp.TM.r_final;            % final drive ratio (dimensionless)

% External Conditions
g = VehData.comp.Ext.g;
rho_air = VehData.comp.Ext.rho_air;

% Driving Cycle: Tq_wheel, Sp_wheel
load('Cycle_short.mat');
del_t = 1;                  % (seconds)
cycle_vel = sch_cycle;
time = [cycle_vel(1,1):del_t:cycle_vel(end,1)];

% Manipulating the driving cycle to get required wheel speed (rad/s) and
% wheel torque (Nm)
v_veh = interp1(cycle_vel(:,1),cycle_vel(:,2),time,'linear')';
cycle_grade = [0 0; time(end) 0];
road_grade = interp1(cycle_grade(:,1),cycle_grade(:,2),time,'linear')';
v_veh_sub = [v_veh; v_veh(end)];
v_veh_sub(1) = [];
v_veh_add = [v_veh(1); v_veh];
v_veh_add(end) = [];
v_acc = (v_veh_sub-v_veh_add)/(2*del_t);

Torque_wheel_temp =   Rtire * (0.5*CD*FA*rho_air*v_veh.*v_veh ...
                                + m_veh*g* cos(road_grade*pi/180).*(Rr1+Rr2.*v_veh+Rr3.*v_veh.*v_veh).*ones(size(v_veh)).*(1+sign(v_veh-0.1))/2 ...
                                + m_veh * v_acc ...
                                + m_veh*g* sin(road_grade*pi/180)).*ones(size(v_veh)).*(1+sign(v_veh-0.1))/2;

Speed_wheel = v_veh/Rtire;
Torque_wheel = Torque_wheel_temp/Gear_eff.*(Torque_wheel_temp>=0) + Torque_wheel_temp*Gear_eff.*(Torque_wheel_temp<0);

OutputTorque = Torque_wheel;       % array of size (T x 1) (Nm) 
OutputSpeed = Speed_wheel;         % array of size (T x 1) (rad/s)

%% 2. Setting time-invariant variables ex) calculating d(m_fc)/dt for points in search space

% define a grid with a set resolution to interpolate the given scattered
% data over the grid (resolution: speed = 1 rad/s, torque = 1 N/m)
grid_sp = round(spE(1)):1:round(spE(end));
grid_tq = tqE(1):1:tqE(end);
[sp_interp, tq_interp] = meshgrid(grid_sp, grid_tq);
sp_interp = sp_interp(:);
tq_interp = tq_interp(:);

% first logical array to check for feasible ENGINE operating points
logicE = (tq_interp' <= repelem(interp1(sp_fullE, tq_maxE, grid_sp),111));

fc_interp = interp2(spE, tqE, fc, sp_interp, tq_interp);
%fc_interp = reshape(fc_interp, [length(grid_tq), length(grid_sp)]);
%mesh(grid_sp, grid_tq, fc_interp);

% V_oc partial gradient w.r.t SOC / R_oc partial gradient w.r.t SOC
grad_V = gradient(Vol_dis, 0.1);
grad_Rcha = gradient(Res_cha, 0.1);
grad_Rdis = gradient(Res_dis, 0.1);

%% 3. Shooting Method

SOC_k = zeros(3, length(OutputSpeed));
p_k = zeros(3, length(OutputSpeed));
Pbatt_k = zeros(3, length(OutputSpeed));
Teng_k = zeros(3, length(OutputSpeed));
Seng_k = zeros(3, length(OutputSpeed));
TMG1_k = zeros(3, length(OutputSpeed));
SMG1_k = zeros(3, length(OutputSpeed));
TMG2_k = zeros(3, length(OutputSpeed));
SMG2_k = zeros(3, length(OutputSpeed));
fc_k = zeros(3, length(OutputSpeed));

SOC_k(:,1) = 0.6;
alpha = 2000;
%alpha = input("please type in value for alpha");

perturb = 0.5;
p_k(1,1) = input('Set an initial guess for the costate at time step 0 (appropriate guess around -300)\n');

p_k(2,1) = p_k(1,1) + perturb;    %perturbed in +ve direction
p_k(3,1) = p_k(1,1) - perturb;    %perturbed in -ve direction

% --------- used for archiving -----------------------
counter = 1;
p0_counter = p_k(1,1);
deltaANDp = fopen('delta_and_p.txt', 'w');
opt_param = [];                                       % to be used as initial setting points for the direct transcription problem
results_alltime = [];                                 % will store SOC_k, Teng_k, Seng_k, SMG1_k, TMG1_k, SMG2_k, TMG2_k (for the entire driving cycle) for the final iteration
results_boundary = [];                                % will store the initial costate, the final costate error, and the total cost at every iteration
% ----------------------------------------------------

while 1
   for t = 1:length(OutputSpeed) - 1
       %% 4. PMP step: finding optimal control
       T_req = OutputTorque(t);
       S_req = OutputSpeed(t);
       
       for lambda = 1:3
           % a) Establishing Search Space: calculating p * f(SOC, u) for points in search space

            SOC_0 = SOC_k(lambda,t);
            p0 = p_k(lambda,t);
            
            %1
            T_MG12 = -1 / (1+SR) * [0, 1; 1+SR, SR] * [repelem(-1 / r_final * T_req, length(tq_interp)); tq_interp']; %first row: T_MG1 for every point in search space
            S_MG12 = [-SR, 1+SR; 1, 0] * [repelem(r_final * S_req, length(sp_interp)); sp_interp'];

            % second/ third logical array to check for feasible MG1/MG2 speeds
            % (speed-wise)
            logicS_MG1 = (S_MG12(1,:) < sp1(end)) & (S_MG12(1,:) > sp1(1));
            logicS_MG2 = (S_MG12(2,:) < sp2(end)) & (S_MG12(2,:) > sp2(1));

            gridT_MG12 = -1 / (1+SR) * [0, 1; 1+SR, SR] * [repelem(-1 / r_final * T_req, length(grid_tq)); grid_tq];
            gridS_MG12 = [-SR, 1+SR; 1, 0] * [repelem(r_final * S_req, length(grid_sp)); grid_sp];

            % fourth / fifth logical array to check for feasible MG1/MG2 torques
            T_MG1max = repelem(interp1(sp_full1, tq_max1, gridS_MG12(1,:)), 111);
            T_MG1min = repelem(interp1(sp_full1, tq_min1, gridS_MG12(1,:)), 111);
            T_MG2max = repelem(interp1(sp_full2, tq_max2, gridS_MG12(2,:)), 111);
            T_MG2min = repelem(interp1(sp_full2, tq_min2, gridS_MG12(2,:)), 111);
            logicT_MG1 = (T_MG12(1,:) < T_MG1max) & (T_MG12(1,:) > T_MG1min);
            logicT_MG2 = (T_MG12(2,:) < T_MG2max) & (T_MG12(2,:) > T_MG2min);

            % calculating required battery power from the electric motor states
            eff1_interp = interp2(sp1, tq1, eff1, S_MG12(1,:), T_MG12(1,:));
            eff2_interp = interp2(sp2, tq2, eff2, S_MG12(2,:), T_MG12(2,:));
            P_batt = eff1_interp .* T_MG12(1,:) .* S_MG12(1,:) + eff2_interp .* T_MG12(2,:) .* S_MG12(2,:);

            % interpolating V_oc at current SOC (V_oc function of only SOC - so we
            % replicate the same values for all points in search space)
            V_oc = repelem(interp1(SOC_ind, Vol_dis, SOC_0), length(sp_interp));
            
            % interpolating R_oc at current SOC (R_oc is function of both SOC and
            % P_batt: sign of P_batt determines which R_oc curve is used to
            % interpolate)
            R_oc = zeros(1, length(sp_interp));
            R_dis = interp1(SOC_ind, Res_dis, SOC_0);
            R_cha = interp1(SOC_ind, Res_cha, SOC_0);

            for i = 1:length(sp_interp)
                if ~isnan(P_batt(i))
                    if P_batt(i) >= 0
                        R_oc(i) = R_dis;
                    else
                        R_oc(i) = R_cha;
                    end
                end
            end
            I_batt = (V_oc - sqrt(V_oc.^2 - 4*P_batt.*R_oc)) ./ (2*R_oc);

            % sixth logical array to check for feasible currents
            I_min = interp1(SOC_ind, Cur_lim_cha, SOC_0);
            I_max = interp1(SOC_ind, Cur_lim_dis, SOC_0);
            logicI = (I_batt < I_max) & (I_batt > I_min);

            f = -I_batt / Q_cap;
            
            %DEBUG ARRAY OF ALL GRIDPOINTS
            %data = [sp_interp'; tq_interp'; fc_interp'; S_MG12; T_MG12; eff1_interp; eff2_interp; logicS_MG1; logicS_MG2; logicT_MG1; logicT_MG2; P_batt; f];
            
            % b) Pontryagin's Minimum Principle (PMP): finding optimal control
            % Viable operating points for consideration for ARGMIN operation is as
            % follows:

            feasible = logicE & logicS_MG1 & logicS_MG2 & logicT_MG1 & logicT_MG2 & logicI;

            H = fc_interp' + p0 * f;
            H_feasible = H(feasible);               % H_feasible: array 1 x 6477 (for this particular case)
            index_feasible = find(feasible);        % index_feasible: array 1 x 6477 of indices of H that are defined
            [H_min, index_min] = min(H_feasible);
            i_opt = index_feasible(index_min);

            Pbatt_k(lambda,t) = P_batt(i_opt);
            Teng_k(lambda,t) = tq_interp(i_opt);
            Seng_k(lambda,t) = sp_interp(i_opt);
            fc_k(lambda,t) = fc_interp(i_opt);
            TMG1_k(lambda,t) = T_MG12(1,i_opt);
            TMG2_k(lambda,t) = T_MG12(2,i_opt);
            SMG1_k(lambda,t) = S_MG12(1,i_opt);
            SMG2_k(lambda,t) = S_MG12(2,i_opt);
       end
       
       fprintf('PMP no. %d: time %d (s) | SOC = %4.4f / %4.4f | p = %7.4f / %7.4f\n', counter, t, SOC_k(2,t), SOC_k(3,t), p_k(2,t), p_k(3,t));
       %% 5. forward integration from step t to step t+1
       V_k = interp1(SOC_ind, Vol_dis, SOC_k(:,t));         %interpolating Voltage at the current time step for perturbed and unperturbed SOCs
       gradV_k = interp1(SOC_ind, grad_V, SOC_k(:,t));      %interpolating Voltage gradient at the current time step for perturbed and unperturbed SOCs
       
       R_kcha = interp1(SOC_ind, Res_cha, SOC_k(:,t));
       R_kdis = interp1(SOC_ind, Res_dis, SOC_k(:,t));
       gradRcha_k = interp1(SOC_ind, grad_Rcha, SOC_k(:,t));
       gradRdis_k = interp1(SOC_ind, grad_Rdis, SOC_k(:,t));
       
       for lambda = 1:3
           if Pbatt_k(lambda,t) >= 0
               R_k = R_kdis(lambda);
               gradR_k = gradRdis_k(lambda);
           else
               R_k = R_kcha(lambda);
               gradR_k = gradRcha_k(lambda);
           end
           %preliminaries
           Num = (sqrt(V_k(lambda)^2 - 4*Pbatt_k(lambda,t)*R_k) - V_k(lambda));
           Denom = (2*Q_cap*R_k);
           gradNum = 1 / (2*sqrt(V_k(lambda)^2 - 4*Pbatt_k(lambda,t)*R_k)) * (2*V_k(lambda)*gradV_k(lambda) - 4*gradR_k*Pbatt_k(lambda,t)) - gradV_k(lambda);
           gradDenom = 2*Q_cap*gradR_k;
           
           %SOC forward integration
           SOC_k(lambda,t+1) = SOC_k(lambda,t) + Num / Denom;
           
           %p forward integration
           part_devf = (gradNum * Denom - Num * gradDenom) / (Denom^2);
           p_k(lambda,t+1) = p_k(lambda,t)*(1 - part_devf);
       end
   end
   total_consum = sum(fc_k(1,:));
   cost = alpha * (SOC_k(1,end) - SOC_k(1,1))^2 + total_consum;
   param_temp = [SOC_k(1,2:end), Teng_k(1,1:end-1), Seng_k(1,1:end-1)];
   opt_param = [opt_param; param_temp];
   
   %% 6. Update Rule
   
   % checking whether integrated p(tf) equals to the one calculated by
   % boundary condition
   
   delta = p_k(1,end) - 2*alpha*(SOC_k(1,end) - SOC_k(1,1));
   results_boundary_temp = [delta, p_k(1,1), cost];
   results_boundary = [results_boundary; results_boundary_temp];
   fprintf("delta: %7.5f\n", delta);
   if (abs(delta) < 0.5)
       %shooting method succeeded in finding correct p0
       results_alltime= [SOC_k(1,2:end); Teng_k(1,1:end-1); Seng_k(1,1:end-1); TMG1_k(1,1:end-1); SMG1_k(1,1:end-1); TMG2_k(1,1:end-1); SMG2_k(1,1:end-1); fc_k(1,1:end-1)];
       fprintf(deltaANDp, '(PMP no.%d) delta: %7.5f | p0_current: %7.5f | total cost: %7.4f\n', counter, delta, p_k(1,1), cost); 
       fclose(deltaANDp);
       break
   else
      %requires update rule to update p0: Newton's method
      %slope = (p_k(3,end) - p_k(2,end)) / (p_k(3,1) - p_k(2,1)); %approximation of slope i.e. derivative of final costate with respect to the initial costate (refer to Donald E. Kirk, 1970, p345)
      Pp = (p_k(3,end) - p_k(2,end)) / (p_k(3,1) - p_k(2,1)); % sensitivity of final costate with respect to perturbations in initial costate
      Px = (SOC_k(3,end) - SOC_k(2,end)) / (p_k(3,1) - p_k(2,1)); % sensitivity of final state with respect to perturbations in initial costate
      slope = 2*alpha*Px - Pp;
      tau = 1;    %learning rate
      p0_new = p_k(1,1) + tau*slope^(-1) * delta;
      del_p0 = tau*slope^(-1) * delta;
      fprintf("slope: %7.5f\n", slope);
      fprintf("del_po: %7.5f\n", del_p0);
      fprintf("p0_new: %7.5f\n", p0_new);
      fprintf('total cost: %7.4f\n', cost);
      fprintf(deltaANDp, '(PMP no.%d) delta: %7.5f | p0_new: %7.5f | slope: %7.5f |total cost: %7.4f\n', counter, delta, p0_new, slope, cost); 
      pause(3)
      p_k(1,1) = p0_new;
      p_k(2,1) = p0_new + perturb;
      p_k(3,1) = p0_new - perturb;
      p0_counter = [p0_counter, p0_new];
      counter = counter + 1;
   end    
end