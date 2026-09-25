%% Rocket Engine Nozzle & Chamber Sizing Script
% Based on MIT Rocket Team Nozzle & Combustion Chamber Design Guide
clear; clc; close all;
%%
%1in diameter
%50psi chamber pressure
%throat will be too small, no reliable manufacturing
%look into combustion chamber general shape and geom, then look into how we
%would design the pintle injectors, look  into thickness so stresses in the
%engine and nozzle
%% 1. Inputs & Operating Conditions
% Propellant properties
OF_ratio  = 2.1;             % O/F Ratio (N2O : Fuel)
w_ethanol = 0.70;            % 70% Ethanol
w_water   = 0.30;            % 30% Water
% Call custom exhaust calculation function (Ensure this is in your path)
[MW, R_spec] = calculateExhaustProps(OF_ratio, w_ethanol, w_water);
% Geometry and Operating Conditions
r_chamber = 0.04286;         % Combustion chamber radius [m]
Ac        = pi * r_chamber^2;% Chamber cross-sectional area [m^2]
Gamma     = 1.2383;          % Specific heat ratio
M_chamber = 0.05;            % Chamber Mach number
temp      = 2151;            % Operating chamber temperature [K]
Lstar     = 0.76;            % Characteristic length L* [m] (Vc / At)
Pc        = 2.758e6;         % Chamber pressure [Pa] (~400 psi)
Pe        = 101352.9;        % Exit pressure [Pa] (~14.7 psi / sea level)
Ru        = 8.314462618;     % Universal gas constant [J/(mol*K)]
g0        = 9.806650;        % Standard acceleration due to gravity [m/s^2]
N_per_lbf = 4.448222;        % Force conversion factor [N/lbf]
%% 2. Reusable Functions
% Isentropic Area Ratio Function (A / A*) based on Mach number and Gamma
calcAreaRatio = @(M, G) ((G + 1)/2)^-((G + 1)/(2*(G - 1))) * ...
                       ((1 + (G - 1)/2 * M^2)^((G + 1)/(2*(G - 1)))) / M;
%% 3. Chamber & Throat Calculations
% Throat Area Calculation
AreaRatio_Chamber = calcAreaRatio(M_chamber, Gamma); % (Ac / A*)
Astar             = Ac / AreaRatio_Chamber;          % Throat area A* [m^2]
% Chamber Volume & Geometry
Vc            = Lstar * Astar;    % Total chamber volume (Vc) [m^3]
ChamberLength = Vc / Ac;          % Cylindrical chamber length [m]
%% 4. Nozzle Exit Calculations
Pratio = Pc / Pe; % Pressure ratio (drives the nozzle expansion)
% Calculate Exit Mach Number algebraically
M_exit = sqrt((2 / (Gamma - 1)) * (Pratio^((Gamma - 1) / Gamma) - 1));
% Calculate sym for fun but commented out
% syms Mach
% assume(Mach>0);
% eqn = Pratio == (1 + ((Gamma - 1)/2) * Mach^2)^(Gamma / (Gamma - 1));
% solM = solve(eqn, Mach);
% Mach = double(solM)
% Exit Area Calculation
AreaRatio_Exit = calcAreaRatio(M_exit, Gamma);       % (Ae / A*)
Ae             = Astar * AreaRatio_Exit;             % Exit area Ae [m^2]
%% 5. Engine Performance Parameters
% Thrust Coefficient (Cf)
Cf = sqrt(((2 * Gamma^2) / (Gamma - 1)) * (2 / (Gamma + 1))^((Gamma + 1) / (Gamma - 1)) * (1 - (Pe / Pc)^((Gamma - 1) / Gamma)));
% Characteristic Velocity (C*)
UpperGamma = sqrt(Gamma * ((Gamma + 1) / 2)^(-(Gamma + 1) / (Gamma - 1)));
Cstar      = (sqrt(R_spec * temp)) / UpperGamma;  % [m/s]
%% 6. Radii Summary
r_throat = sqrt(Astar / pi); % Throat radius [m]
r_exit   = sqrt(Ae / pi);    % Exit radius [m]
%% 7. Results Display
fprintf('\n=== GEOMETRY & SIZING ===\n');
fprintf('Chamber Radius (Rc) : %.6f m\n', r_chamber);
fprintf('Throat Radius (Rt)  : %.6f m\n', r_throat);
fprintf('Exit Radius (Re)    : %.6f m\n', r_exit);
fprintf('---------------------------------\n');
fprintf('Chamber Area (Ac)   : %.6f m^2\n', Ac);
fprintf('Throat Area (A*)    : %.6f m^2\n', Astar);
fprintf('Exit Area (Ae)      : %.6f m^2\n', Ae);
fprintf('Chamber Volume (Vc) : %.6f m^3\n', Vc);
fprintf('Chamber Length      : %.6f m\n', ChamberLength);
fprintf('Area Ratio (Ac/A*)  : %.6f\n', AreaRatio_Chamber);
fprintf('Expansion Ratio (Ae/A*): %.6f\n', AreaRatio_Exit);
fprintf('\n=== PERFORMANCE & FLOW ===\n');
fprintf('Pressure Ratio (Pc/Pe): %.6f\n', Pratio);
fprintf('Exit Mach Number (Me) : %.6f\n', M_exit);
fprintf('Thrust Coeff (Cf)     : %.6f\n', Cf);
fprintf('Upper Gamma           : %.6f\n', UpperGamma);
fprintf('Char. Velocity (C*)   : %.6f m/s\n', Cstar);
%% figure out length of each section and shape of the curve
% all isentropic are reversible process, mirror nozzle, irl it isn't
% isentropic
% on the wall there are friction and heat flux and boundary layers, so lose
% some energy, but basically isentropic
%The length of the nozzle isn't taken into account, so look into length,
% if too short it starts being not quasi-1D, avoid shocks, figure out how
% to mesh with combustion chamber, want to be fully combusted before the
% nozzle, few equations on the MIT
% most basic nozzle is straight line, throat too sharp causes a lot of
% heat, we want a parabolic change for a sharper transition at the
% beginning and end of the nozzle, but a gradual change at the middle of
% the nozzle
%if we can get a gerneral mass flow and thrust + basic length, can start
%doing cfd and see what we get
%% 8. MIT Geometry Steps: Converging & Diverging Sections
% Solves for the lengths mentioned in the CFD notes above
% --- Converging Section & Final Cylinder Cut Length ---
theta_deg  = 30;             % Contraction half-angle [degrees]
R_arc_conv = 1.5 * r_throat; % Throat entrance arc radius
% Calculate Converging Length (L_conv)
L_conv = (r_chamber - r_throat + R_arc_conv * (secd(theta_deg) - 1)) / tand(theta_deg);
% Calculate Converging Volume (Truncated Cone Formula)
V_conv = (1/3) * pi * L_conv * (r_chamber^2 + r_chamber * r_throat + r_throat^2);
% Extract True Cylindrical Chamber Length
V_cyl = Vc - V_conv;         % Remaining volume for straight pipe
L_cyl = V_cyl / Ac;          % True Cylinder cut length [m]
% --- Diverging Section Length ---
alpha_deg  = 15;             % Divergent half-angle [degrees]
R_arc_div  = 1.5 * r_throat; % Throat exit arc radius
% Calculate Diverging Nozzle Length (L_div)
L_div = (r_exit - r_throat + R_arc_div * (secd(alpha_deg) - 1)) / tand(alpha_deg);
%% 9. Manufacturing Lengths Display
fprintf('\n=== MANUFACTURING LENGTHS & VOLUMES ===\n');
fprintf('Total Chamber Vol (Vc) : %.6f m^3 (L* = %.2f m)\n', Vc, Lstar);
fprintf('Idealized Tank Length  : %.6f m   (From basic Vc/Ac)\n', ChamberLength);
fprintf('---------------------------------------------------\n');
fprintf('1. True Cylinder Length: %.6f m  | Vol: %.6f m^3\n', L_cyl, V_cyl);
fprintf('2. Converging Length   : %.6f m  | Vol: %.6f m^3 (theta = %d deg)\n', L_conv, V_conv, theta_deg);
fprintf('3. Diverging Length    : %.6f m  |                 (alpha = %d deg)\n', L_div, alpha_deg);
fprintf('---------------------------------------------------\n');
fprintf('Total Engine Length    : %.6f m  (Injector face to Nozzle exit)\n\n', L_cyl + L_conv + L_div);
%% Calling the geom
% Call nozzle-only contour:
generateNozzleContour(r_throat, r_exit);
% Or call full engine contour:
generateFullEngineContour(r_chamber, r_throat, r_exit, L_cyl);
% Specify your desired wall thickness in mm (e.g., 3.0 mm)
wall_thickness_mm = 3.0;
generateParabolicEngineContour(r_chamber, r_throat, r_exit, L_cyl, wall_thickness_mm);
%% Finding mass flow rate, assuming 100% C* efficiency
% Experimental used in MIT calcs - C* = (At * Pc) / m*
% Real C* would be 88-95% efficiency on average
% Cstar_real = Cstar * .9
mtotal = (Astar * Pc) / Cstar;
% oxidizer and fuel amounts based on our current 2.1:1 O/F ratio
mfuel = mtotal / (1 + OF_ratio);
mox = mtotal - mfuel;
% results display
fprintf('---------------------------------------------------\n')
fprintf('=== MASS FLOW RATES ===\n')
fprintf('Total Mass Flow: %.6f kg/s\n', mtotal);
fprintf('Fuel Mass Flow:  %.6f kg/s\n', mfuel);
fprintf('Oxidizer Mass Flow: %.6f kg/s\n', mox);
%% Isp and Thrust Calculations
% These use the equations given from the MIT site but algebraically change
% them to fit our needs
% Specific Impulse Calculation
% Evaluated directly from combustion efficiency (c*) and nozzle expansion (Cf)
Isp            = (Cf * Cstar) / g0;                 % Specific Impulse [s]
c_eff          = Isp * g0;                          % Effective exhaust velocity [m/s]
% Total Thrust Calculations
% Method A: Calculated via Mass Flow Rate & Effective Exhaust Velocity
F_total_N      = mtotal * c_eff;                  % Total Thrust [N]
% Method B: Direct verification via Pressure & Area (F = Cf * Pc * At)
% F_total_N    = Cf * P_c * A_t;                    
% Convert Thrust to Imperial Units
F_total_lbf    = F_total_N / N_per_lbf;             % Total Thrust [lbf]
fprintf('\n=== ENGINE PERFORMANCE OUTPUT ===\n');
fprintf('Specific Impulse (Isp):      %8.6f s\n',    Isp);
fprintf('Effective Exhaust Velocity:  %8.6f m/s\n',  c_eff);
fprintf('Total Thrust (N):            %8.6f N\n',    F_total_N);
fprintf('Total Thrust (lbf):          %8.6f lbf\n',  F_total_lbf);
fprintf('---------------------------------------------------\n\n');
