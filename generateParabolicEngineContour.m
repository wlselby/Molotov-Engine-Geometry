function generateParabolicEngineContour(r_chamber, r_throat, r_exit, L_cyl, wall_thickness)
   % GENERATECONTOUREDIDANDOD Generates ID and OD profiles sharing 100%
   % identical X-coordinates so both profiles end at the exact same exit plane.
   %
   % Inputs:
   %   r_chamber      : Chamber radius [m]
   %   r_throat       : Throat radius [m]
   %   r_exit         : Nozzle exit radius [m]
   %   L_cyl          : Chamber cylinder length [m]
   %   wall_thickness : Wall thickness in millimeters [mm] (e.g., 3.0)
   if nargin < 5
       wall_thickness = 3.0; % Default wall thickness in mm
   end
   %% 1. Design Parameters
   theta_conv   = 30;                 % Converging half-angle [deg]
   R_arc_cham   = 1.5 * r_throat;     % Smooth chamber entrance arc radius [m]
   R_arc_conv   = 1.5 * r_throat;     % Throat entrance arc radius [m]
   R_arc_div    = 0.382 * r_throat;   % 80% Bell upstream arc radius [m]
   theta_n      = 20.0;                % Initial parabola angle [deg]
   theta_e      = 7.0;                 % Final exit angle [deg]
   %% 2. Calculate Key Coordinates (Smooth Tangencies)
   % Throat Upstream Tangency
   x_conv_tan = -R_arc_conv * sind(theta_conv);
   y_conv_tan = (r_throat + R_arc_conv) - R_arc_conv * cosd(theta_conv);
   % Chamber Entrance Arc Tangency (Blends Cylinder to 30-deg Cone)
   y_cham_tan   = (r_chamber - R_arc_cham) + R_arc_cham * cosd(theta_conv);
   x_cham_tan   = x_conv_tan - (y_cham_tan - y_conv_tan) / tand(theta_conv);
   x_cham_start = x_cham_tan - R_arc_cham * sind(theta_conv);
  
   x_injector   = x_cham_start - L_cyl;
   % Diverging Bell Coordinates
   x_N = R_arc_div * sind(theta_n);
   y_N = (r_throat + R_arc_div) - R_arc_div * cosd(theta_n);
   L_cone = (r_exit - r_throat) / tand(15);
   L_bell = 0.80 * L_cone;
   x_E = L_bell;
   y_E = r_exit;
   m1 = tand(theta_n);
   m2 = tand(theta_e);
   x_Q = (y_E - y_N + m1*x_N - m2*x_E) / (m1 - m2);
   y_Q = y_N + m1 * (x_Q - x_N);
   %% 3. Generate Inner Gas Path Segments (ID)
   num_pts = 60;
  
   % Segment 0: Straight Chamber Cylinder
   x0 = linspace(x_injector, x_cham_start, num_pts);
   y0 = ones(1, num_pts) * r_chamber;
   % Segment 1: Chamber Entrance Arc (0 -> -30 deg)
   angles_cham = linspace(0, deg2rad(theta_conv), num_pts);
   x1 = x_cham_start + R_arc_cham * sin(angles_cham);
   y1 = r_chamber - R_arc_cham * (1 - cos(angles_cham));
   % Segment 2: Converging Cone (-30 deg)
   x2 = linspace(x_cham_tan, x_conv_tan, num_pts);
   y2 = linspace(y_cham_tan, y_conv_tan, num_pts);
   % Segment 3: Throat Entrance Arc (-30 -> 0 deg)
   angles_conv = linspace(-deg2rad(theta_conv), 0, num_pts);
   x3 = R_arc_conv * sin(angles_conv);
   y3 = (r_throat + R_arc_conv) - R_arc_conv * cos(angles_conv);
   % Segment 4: Throat Exit Arc (0 -> +20 deg)
   angles_div = linspace(0, deg2rad(theta_n), num_pts);
   x4 = R_arc_div * sin(angles_div);
   y4 = (r_throat + R_arc_div) - R_arc_div * cos(angles_div);
   % Segment 5: Parabolic Bell Curve (Quadratic Bézier)
   t  = linspace(0, 1, num_pts);
   x5 = (1 - t).^2 * x_N + 2*(1 - t).*t * x_Q + t.^2 * x_E;
   y5 = (1 - t).^2 * y_N + 2*(1 - t).*t * y_Q + t.^2 * y_E;
   % Combine Segments (Meters to Millimeters)
   X_in_m = [x0, x1(2:end), x2(2:end), x3(2:end), x4(2:end), x5(2:end)]';
   Y_in_m = [y0, y1(2:end), y2(2:end), y3(2:end), y4(2:end), y5(2:end)]';
  
   X_mm    = X_in_m * 1000;
   Y_in_mm = Y_in_m * 1000;
   %% 4. Calculate OD with Identical X-Coordinates
   dX = gradient(X_mm);
   dY = gradient(Y_in_mm);
   theta_local = atan2(dY, dX); % Local wall angle along the profile
   X_out_mm = X_mm; % Identical X length!
   Y_out_mm = Y_in_mm + (wall_thickness ./ cos(theta_local)); % Correct normal thickness
   %% 5. Export Files
   % --- File 1: Inner Profile (ID) ---
   matrix_ID = [X_mm, Y_in_mm, zeros(length(X_mm), 1)];
   writematrix(matrix_ID, 'Engine_Inner_Profile_ID.txt', 'Delimiter', 'tab');
   % --- File 2: Outer Profile (OD) ---
   matrix_OD = [X_out_mm, Y_out_mm, zeros(length(X_out_mm), 1)];
   writematrix(matrix_OD, 'Engine_Outer_Profile_OD.txt', 'Delimiter', 'tab');
   %% 6. Console Status Report
   fprintf('\n=== FILES EXPORTED (MATCHING LENGTHS) ===\n');
   fprintf('Injector Plane X : %.2f mm (Both ID and OD)\n', X_mm(1));
   fprintf('Nozzle Exit Plane X: %.2f mm (Both ID and OD)\n', X_mm(end));
   fprintf('Wall Thickness   : %.2f mm\n\n', wall_thickness);
end

