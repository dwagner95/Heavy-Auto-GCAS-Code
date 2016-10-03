% c_130J

% \module genesis

% Initial and Reference Conditions
GS_init_v_fps = 295.0; %Initial true airspeed (ft/sec)
GS_alpha_ref = 7.5; %Trim gain for winds (nondim)
GS_Vt_ref = 295.0; %True airspeed where aero coefficients are valid (ft/sec)

% Aircraft Dimensions
GS_wing_area = 6200.0; %Wing area (ft^2)
GS_cbar = 30.93; %Mean aerodynamic chord (ft)
GS_b = 219.17; %Wing span (ft)
GS_pilot_x = 81.0; %Pilot station X-axis location (ft fwd of CG)
GS_pilot_y = 0.0; %Pilot station Y-axis location (ft right of CG)
GS_pilot_z = -25.0; %Pilot station Z-axis location (ft below CG)

% Mass Properties
GS_cg_ref = 0.25; %CG where Cmalpha is valid (frac. MAC)
GS_weight = 728000.0; %Gross weight (lb)
GS_Ixx = 32500000.0; %Moment of inertia about X-axis (slug-ft^2)
GS_Iyy = 32800000.0; %Moment of inertia about Y-axis (slup-ft^2)
GS_Izz = 61500000.0; %Moment of inertia about Z-axis (slup-ft^2)
GS_Ixz = -3340000.0; %Product of inertia for X and Z axes (slup-ft^2)

% Lift Coefficients
CL0 = 0.25; %CL at alpha = 0 (nondim)
CLalpha = 5.37; %CL per alpha (rad^-1)
CLq = 0.0; %CL per pitch rate ((rad/sec)^-1) 
CLalpha_dot = 0.0; %CL per alpha rate ((rad/sec)^-1)
CLdelta_e = 0.324; %CL per elevator deflection (rad^-1)
CLdelta_s = 0.324; %CL per stabilizer deflection (rad^-1)
CLdelta_f = 0.5; %CL per flap deflection (rad^-1)
CLu = 0.0; %CL per delta velocity ((ft/sec)^-1)
CLmax = 10000.0; %Maximum CL (nondim)

% Pitching Moment Coefficients
Cm0 = 0.0; %Cm for fuselage (nondim)
Cmalpha = -1.57; %Cm per alpha (rad^-1)
Cmq = -24.1; %Cm per pitch rate ((rad/sec)^-1)
Cmalpha_dot = -7.8; %Cm per alpha rate ((rad/sec)^-1)
Cmdelta_e = -1.35; %Cm per elevator deflection (rad^-1)
Cmdelta_s = -1.35; %Cm per stabilizer deflection (rad^-1)
Cmdelta_f = 0.0; %Cm per flap deflection (rad^-1)
Cmdelta_sb = 0.0; %Cm per speedbrake deflection (rad^-1)
Cmdelta_g = 0.0; %Cm per gear deflection (nondim)
Cmu = 0.0; %Cm per delta velocity ((ft/sec)^-1)

% Drag Coefficients
CD0 = 0.025; %CD at CL = 0 (nondim)
CDq = 0.0; %CD per pitch rate ((rad/sec)^-1)
CDalpha = 0.0; %Slope of parabolic drag polar (nondim)
CDCL2 = 0.018; %CD per alpha (rad^-1)
CDdelta_e = 0.14; %CD per elevator deflection (rad^-1)
CDdelta_s = 0.14; %CD per stabilizer deflection (rad^-1)
CDdelta_f = 0.0; %CD per flap deflection (rad^-1)
CDdelta_g = 0.01; %CD per gear deflection (nondim)
CDdelta_sb = 0.01; %CD per speedbrake deflection (rad^-1)
CDu = 0.0; %CD per delta velocity ((ft/sec)^-1)

% Sideforce Coefficients
CY0 = 0.0; %CY at beta = 0 (non-dimensional)
CYp = 0.0; %CY per roll rate ((rad/sec)^-1)
CYr = 0.0; %CY per yaw rate ((rad/sec)^-1)
CYbeta = -0.742; %CY per beta (rad^-1)
CYdelta_a = 0.0; %CY per aileron deflection (rad^-1)
CYdelta_r = 0.2065; %CY per rudder deflection (rad^-1)

% Rolling Moment Coefficients
Clp = -0.41; %Cl per roll rate ((rad/sec)^-1)
Clr = 0.309; %Cl per yaw rate ((rad/sec)^-1)
Clbeta = -0.14; %Cl per beta (rad^-1)
Cldelta_a = -0.053; %Cl per aileron deflection (rad^-1)
Cldelta_r = -0.0157; %Cl per rudder deflection (rad^-1)

% Yawing Moment Coefficients
Cnp = -0.1110; %Cn per roll rate ((rad/sec)^-1)
Cnr = -0.225; %Cn per yaw rate ((rad/sec)^-1)
Cnbeta = 0.149; %Cn per beta (rad^-1)
Cndelta_a = -0.0032; %Cn per aileron deflection (rad^-1)
Cndelta_r = -0.1060; %Cn per rudder deflection (rad^-1)

% Landing Gear Characteristics
gear_main_x = -10.0; %Main gear X-axis location (ft fwd of CG)
gear_main_y = 18.0; %Right main gear Y-axis location (ft right of CG)
gear_main_z = 25.5; %Main gear Z-axis location (ft below CG)
gear_nose_x = 70.0; %Nose gear X-axis location (ft fwd of CG)
gear_nose_z = 25.0; %Nose gear Z-axis location (ft below CG)
gear_main_spring_constant = 600000.0; %Main gear spring constant (lb/ft/sec)
gear_nose_spring_constant = 580000.0; %Nose gear spring constant (lb/ft/sec)
gear_main_damping_constant = 10000.0; %Main gear damper constant (lb/ft/sec)
gear_nose_damping_constant = 2000.0; %Nose gear damper constant (lb/ft/sec)
gear_main_rolling_friction = 0.0012; %Main gear rolling coefficient of friction (non-dimensional)
gear_nose_rolling_friction = 0.0006; %Nose gear rolling coefficient of friction (non-dimensional)

% Engine Characteristics
dx_omega = 0.5; %Engine first-order frequency (1/tau) (rad/sec)
dx_min_power = 750.0; %Minimum (idle) thrust per engine (lb)
dx_max_power = 82200.0; %Maximum thrust per engine (lb)
fuel_flow_max = 4000.0; %Fuel flow at 100% throttle (lb/hr)
dx_delay = 0.0; %Engine time delay (sec)
dx_y = 50.0; %Engine location right of CG (ft)
dx_z = 0.7; %Engine location below CG (ft)
dx_incidence = 0.0122; %Engine incidence above X-Y plane (rad)

% Pitch FCS
dec_per_des = -30.0; %Elevator command per longitudinal stick deflection (deg/%)
dec_per_q = 0.15; %Elevator command per pitch rate (deg/deg/sec)
q_wash_freq = 1.0; %Pitch rate washout filter frequency (rad/sec)
dec_per_alpha = 0.0; %Elevator command per alpha (deg/deg)
dec_per_nz = 0.0; %Elevator command per Nz (deg/G)
pitch_trim_rate = -2.0; %Pitch trim rate (deg/sec)

% Stabilizer
stab_to_trim_pos = -12.8; %Stabilizer takeoff trim position (deg)
ds_trim_max = 30.0; %Stabilizer actuator negative trim limit (deg)
ds_trim_min = -30.0; %Stabilizer actuator positive trim limit (deg)

% Elevator
de_omega = 30.0; %Elevator actuator natural frequency (rad/sec)
de_zeta = 0.7; %Elevator actuator damping ratio (nondim)
de_min_rate = -100.0; %Elevator actuator maximum negative-going rate (deg/sec)
de_max_rate = 100.0; %Elevator actuator maximum positive-going rate (deg/sec)
de_min_position = -30.0; %Elevator actuator negative position limit (deg)
de_max_position = 30.0; %Elevator actuator positive position limit (deg)
de_delay = 0.020; %Elevator actuator time delay (sec)

% Roll FCS
dac_per_das = -30.0; %Aileron command per lateral stick deflection (deg/%)
roll_trim_rate = -2.0; %Roll trim rate (deg/sec)
dac_per_p = 0.265; %Aileron command per roll rate (deg/(deg/sec))

% Aileron
da_omega = 30.0; %Aileron actuator naturla frequency (rad/sec)
da_zeta = 0.7; %Aileron actuator damping ratio (nondim)
da_min_rate = -150.0; %Aileron actuator maximum negative-going rate (deg/sec)
da_max_rate = 150.0; %Aileron actuator maximum positive-going rate (deg/sec)
da_min_position = -30.0; %Aileron actuator negative position limit (deg)
da_max_position = 30.0; %Aileron actuator positive position limit (deg)
da_trim_min = -6.0; %Aileron actuator negative trim authority limit (deg)
da_trim_max = 6.0; %Aileron actuator positive trim authority limit (deg)
da_delay = 0.010; %Aileron actuator time delay (sec)

% Yaw FCS
drc_per_drp = -20.0; %Rudder command per rudder pedal deflection (deg/%)
drc_per_betadot = 0.0; %Rudder command per beta rate (deg/deg/sec)
drc_per_r = 0.15; %Rudder command per yaw rate (deg/deg/sec)
r_wash_freq = 2.0; %Yaw rate washout filter frequency (rad/sec)
drc_per_ay = 0.96; %Rudder command per lateral acceleration (deg/G)
drc_per_dac = 0.69; %Rudder command per aileron command (deg/deg)
yaw_trim_rate = -2.0; %Yaw trim rate (deg/sec)

% Rudder
dr_omega = 30.0; %Rudder actuator natural frequency (rad/sec)
dr_zeta = 0.7; %Rudder actuator damping ratio (nondim)
dr_min_rate = -100.0; %Rudder actuator maximum negative-going rate (deg/sec)
dr_max_rate = 100.0; %Rudder actuator maximum positive-going rate (deg/sec)
dr_min_position = -20.0; %Rudder actuator negative position limit (deg)
dr_max_position = 20.0; %Rudder actuator positive position limit (deg)
dr_trim_min = -20.0; %Rudder actuator negative trim authority limit (deg)
dr_trim_max = 20.0; %Rudder actuator positive trim authority limit (deg)
dr_delay = 0.0; %Rudder actuator time delay

% Secondary Surfaces
dg_extend_time = 2.5; %Gear extension time (sec)
dg_retract_time = 2.5; %Gear retraction time (sec)
dfc_max = 40.0; %Maximum flap command (deg)
df_rate = 8.0; %Flap actuator rate (deg/sec)
dsbc_max = 20.0; %Speedbrake maximum position (deg)
dsbc_ext_rate = 3.0; %Speedbrake extension rate (deg/sec)
dsb_rate = 3.0; %Speedbrake actuator rate (deg/sec)

% Display Parameters
disp_vselected_kts = 170.0; %Selected airspeed (knots)
disp_vref_kts = 170.0; %Reference approach airspeed display (knots)
disp_vmin_kts = 130.0; %Minimum airspeed display (knots)
disp_vmax_kts = 750.0; %Maximum airspeed display (knots)
disp_alpha_max = 26.0; %Position of yellow "angel wings" on B777 PFD (deg)
