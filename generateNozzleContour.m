function [X, Y] = generateNozzleContour(r_throat, r_exit)
   % GENERATENOZZLECONTOUR Generates nozzle contour and saves a .txt file for SolidWorks
  
   %% Design Parameters (MIT Baseline)
   theta_conv = 30;             % Converging section half-angle [deg]
   alpha_div  = 15;             % Diverging section half-angle [deg]
   R_arc_conv = 1.5 * r_throat; % Arc smoothing into the throat
   R_arc_div  = 1.5 * r_throat; % Arc smoothing out of the throat
   %% Calculate Tangency Points
   x_conv_tan = -R_arc_conv * sind(theta_conv);
   y_conv_tan = (r_throat + R_arc_conv) - R_arc_conv * cosd(theta_conv);
   x_div_tan = R_arc_div * sind(alpha_div);
   y_div_tan = (r_throat + R_arc_div) - R_arc_div * cosd(alpha_div);
   x_exit = x_div_tan + (r_exit - y_div_tan) / tand(alpha_div);
   y_exit = r_exit;
   %% Generate Curve Segments
   num_points = 50;
   % Segment 1: Converging Arc
   angles_conv = linspace(-deg2rad(theta_conv), 0, num_points);
   x1 = R_arc_conv * sin(angles_conv);
   y1 = (r_throat + R_arc_conv) - R_arc_conv * cos(angles_conv);
   % Segment 2: Diverging Arc
   angles_div = linspace(0, deg2rad(alpha_div), num_points);
   x2 = R_arc_div * sin(angles_div);
   y2 = (r_throat + R_arc_div) - R_arc_div * cos(angles_div);
   % Segment 3: Diverging Straight Line
   x3 = linspace(x_div_tan, x_exit, num_points);
   y3 = linspace(y_div_tan, y_exit, num_points);
   %% Combine Points into XYZ Matrix
   X = [x1, x2(2:end), x3(2:end)]';
   Y = [y1, y2(2:end), y3(2:end)]';
   Z = zeros(length(X), 1);
   point_matrix = [X, Y, Z];
  
   %% Save as tab-delimited .txt file for SolidWorks
   writematrix(point_matrix, 'Nozzle_Only_Contour.txt', 'Delimiter', 'tab');
   %% Plot Geometry
   figure('Name', 'Nozzle Only Contour (MIT Conical)', 'Color', 'w');
   hold on; grid on; axis equal;
   plot([x_conv_tan - 0.01, x_exit + 0.01], [0, 0], 'k-.', 'LineWidth', 1.5);
   plot(X, Y, 'b-', 'LineWidth', 2);
   plot(X, -Y, 'b-', 'LineWidth', 2);
   xline(0, 'r--', 'Throat', 'LabelVerticalAlignment', 'bottom');
   title('Nozzle Only Fluid Volume (MIT Baseline)');
   xlabel('Length (m)'); ylabel('Radius (m)');
   legend('Centerline', 'Wall Contour', 'Location', 'best');
   fprintf('Success: "Nozzle_Only_Contour.txt" saved to folder.\n');
end

