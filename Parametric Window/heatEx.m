function heatEx
% Define constants
k = 0.6; % Thermal conductivity of water (W/mK)
Cp = 4187; % Specific heat capacity of water (J/kgK)
rho = 1000; % Density of water (kg/m^3)

% Define design variables and ranges
Vc_range = linspace(0, 1000/60, 20); % Flow rate range for cold fluid (m^3/s)
WCH_range = linspace(100e-6, 1000e-6, 20); % Width range for channels (m)
DCH_range = linspace(100e-6, 1000e-6, 20); % Depth range for channels (m)
SCH_range = linspace(100e-6, 500e-6, 20); % Spacing range for channels (m)
LCH_range = linspace(0, 20, 20); % Length range for channels (m)

% Create arrays to store results
Q_array = zeros(length(Vc_range), length(WCH_range), length(DCH_range), length(SCH_range), length(LCH_range));
P_array = zeros(length(Vc_range), length(WCH_range), length(DCH_range), length(SCH_range), length(LCH_range));

% Loop through all combinations of design variables
for i = 1:length(Vc_range)
    for j = 1:length(WCH_range)
        for k = 1:length(DCH_range)
            for l = 1:length(SCH_range)
                for m = 1:length(LCH_range)
                    % Calculate Reynolds number
                    Re = Vc_range(i)*DCH_range(k)/(WCH_range(j)*rho)
                    
                    % Calculate friction factor
                    f = 64/Re;
                    
                    % Calculate heat transfer coefficient
                    h = 2*k/(DCH_range(k)*f);
                    
                    % Calculate pressure drop
                    P = f*LCH_range(m)*rho*Vc_range(i)^2/(2*DCH_range(k));
                    
                    % Calculate heat transfer rate
                    Q = h*WCH_range(j)*DCH_range(k)*LCH_range(m)*(80-40);
                    
                    % Store results in arrays
                    Q_array(i,j,k,l,m) = Q;
                    P_array(i,j,k,l,m) = P;
                end
            end
        end
    end
end

% Find the combination of design variables that satisfies the requirements and constraints
% (i.e. Q >= 0.2 W and P <= 20 Pa)

% Check if there are any valid combinations of design variables
Q_target = 0.2; % Target heat transfer rate (W)
P_target = 20; % Target pressure drop (Pa)

Q_diff = abs(Q_array - Q_target);
P_diff = abs(P_array - P_target);
[Q_idx, P_idx] = find(Q_diff <= 0.01*Q_target & P_diff <= P_target);

if ~isempty(Q_idx)
    % Select the optimal design variables that give the minimum pressure drop for a given heat transfer rate 
    [~, idx] = min(P_array(Q_idx(1), :, P_idx(1), :, :), [], 'all', 'linear');
    [Vc_idx, WCH_idx, DCH_idx, SCH_idx, LCH_idx] = ind2sub(size(P_array), idx);
    
% Print the combination of design variables and corresponding results
 disp(['Vc: ' num2str(Vc_range(Vc_idx)) ' m^3/s']); disp(['WCH: ' num2str(WCH_range(WCH_idx)) ' m']); disp(['DCH: ' num2str(DCH_range(DCH_idx)) ' m']); disp(['SCH: ' num2str(SCH_range(SCH_idx)) ' m']); disp(['LCH: ' num2str(LCH_range(LCH_idx)) ' m']); disp(['Q: ' num2str(Q_array(Vc_idx, WCH_idx, DCH_idx, SCH_idx, LCH_idx)) ' W']); disp(['P: ' num2str(P_array(Vc_idx, WCH_idx, DCH_idx, SCH_idx, LCH_idx)) ' Pa']);
% Plot heat transfer rate as a function of pressure drop for the selected design variables
Q_vs_P = squeeze(Q_array(Vc_idx, WCH_idx, DCH_idx, SCH_idx, :));
P_vs_Q = squeeze(P_array(Vc_idx, WCH_idx, DCH_idx, SCH_idx, :));
figure;
subplot(2,1,1);
plot(P_vs_Q, LCH_range);
xlabel('Pressure drop (Pa)');
ylabel('Length (m)');
title('Length vs Pressure Drop');
subplot(2,1,2);
plot(Q_vs_P, LCH_range);
xlabel('Heat transfer rate (W)');
ylabel('Length (m)');
title('Length vs Heat Transfer Rate');
sgtitle(['Design Variables: Vc = ' num2str(Vc_range(Vc_idx)) ' m^3/s, WCH = ' num2str(WCH_range(WCH_idx)) ' m, DCH = ' num2str(DCH_range(DCH_idx)) ' m, SCH = ' num2str(SCH_range(SCH_idx)) ' m, LCH = ' num2str(LCH_range(LCH_idx)) ' m']);

% Plot pressure drop as a function of length for the selected design variables
P_vs_L = squeeze(P_array(Vc_idx, WCH_idx, DCH_idx, SCH_idx, :));
figure;
plot(LCH_range, P_vs_L);
xlabel('Length (m)');
ylabel('Pressure drop (Pa)');
title('Length vs Pressure Drop');
sgtitle(['Design Variables: Vc = ' num2str(Vc_range(Vc_idx)) ' m^3/s, WCH = ' num2str(WCH_range(WCH_idx)) ' m, DCH = ' num2str(DCH_range(DCH_idx)) ' m, SCH = ' num2str(SCH_range(SCH_idx)) ' m, LCH = ' num2str(LCH_range(LCH_idx)) ' m']);

end
end

