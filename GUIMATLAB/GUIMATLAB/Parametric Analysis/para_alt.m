function [net_thrust_array , sfc_array] = para_alt(h1, h2, mach_1 , mach_2)
num = ((h2 - h1)/1000)+1;
Temp = zeros(6,1);
Press = zeros(6,1);
net_thrust_array = zeros(1,num);
spec_net_thrust_array = zeros(1,num);
sfc_array = zeros(1,num);
comp_ratio_array = zeros(1,num);
air_mass_array = zeros(1,num);
i = 1;
j = 1;

cp_gas = 1000;  %J/Kg.Kelvin
Qhv = 43260000; %J/kg
gamma_gas = 1.4;
R = 287;
h = 10000;
Eff_c = 1;
BET=1050;
PR_burn = 0.96;
Eff_Noz = 1;
A_noz = 0.3;
Eff_isen_turb = .9;
Eff_mechn_turb = 1;


i = 1;
j = 1;

comp_ratio = 4;


for mach = mach_1:.1:mach_2
j = 1;
for height = h1:1000:h2
    Mo = mach
    h =height
    if h <=11000
        Tamb = (15.04-(0.00649*h));
        Pamb = (101.29* ((Tamb+273.1)/288.08)^5.256);
    elseif h >11000 && h <25000
            Tamb = -56.46;
            Pamb = 22.65 * (exp(1.73-(.000157*h)));
        else 
            Tamb = -131.21 + (.00299 * h);
            Pamb = 2.488 * (((Tamb+273.1)/(216.6))^-11.388);
    end
    Velo_sound = (gamma_gas*R*(Tamb + 273))^0.5;
Velo_1 = Mo*Velo_sound;
Dynamic_Temp = (Velo_1^2)/(2*cp_gas);



Pt = Pamb*1000*(1 + Dynamic_Temp/(Tamb + 273))^(gamma_gas/(gamma_gas-1)); %Pascals
Tt = Tamb + 273.16 + Dynamic_Temp; 
%Kelvin

Po = Pt ;%Pascals
To = Tt ;%Kelvin

%Po = Pamb*1000 %Pascals
%To = Tamb + 273.16 %Kelvin
so  = 1391; % This needs to be automated

Temp(1,1) = To;
Press(1,1) = Po;
%s_entro(1,1) = so;
vo=(R*To)/Po;
s_vol(1,1) = vo;
h = height;
    Tamb = (15.04-(0.00649*h))
    Pamb = (101.29* ((Tamb+273.1)/288.08)^5.256)

PR_c = comp_ratio;
height_array(i,j) = height;
Velo_1 = Mo*Velo_sound;

%Compressor Function Call
Temp(2,1) = To;
Press(2,1) = Po;
[ T2 , p2  W_c ] = Comp_Func (  To , Po , PR_c , Eff_c , gamma_gas , cp_gas);
Temp(3,1) = T2;
Press(3,1) = p2;

%Burner Function Call
[ T3 ,p3 ,AF ] = burner_func (  T2 , p2 ,  cp_gas ,PR_burn,BET,Qhv );

Temp(4,1) = T3;
Press(4,1) = p3;

%Turbine Function Call

W_t = W_c;

[ T4 , p4 ] = Turb_Func( T3 , p3 , Eff_isen_turb , Eff_mechn_turb , W_t , gamma_gas , cp_gas);

Temp(5,1) = T4;
Press(5,1) = p4;

%Nozzle Function Call

[ m5 , T5 , p5 , Velo_5 ] = Nozzle_Func( T4 , p4 , Eff_Noz , gamma_gas , A_noz,R);

Temp(6,1) = T5;
Press(6,1) = p5;

m_fuel = m5 /(1+AF);
m_inlet = m5-m_fuel;

%Thrust Function Call
[net_thrust , net_thrust_spec] = thrust_func(Velo_1,Velo_5, m5,Po,p5,A_noz,m_inlet);
net_thrust 
net_thrust_array(i,j) = net_thrust;
spec_net_thrust_array(i,j) = net_thrust_spec;

%SFC Calculation
m_fuel = m5 /(1+AF);
m_inlet = m5-m_fuel
Theeta= To/(273.16+15);
Pressure_Theeta = Po/(101325);
m_sea = m_inlet * (((Theeta)^0.5)/(Pressure_Theeta))
air_mass_array(i,j) = m_inlet;
specific_fuel = m_fuel/net_thrust
sfc_array(i,j) = specific_fuel;
fprintf("-----------------------------------------------------");


j = j + 1;

end



%PLOTS

figure(1)
plot(height_array(i,:), net_thrust_array(i,:) ,'-o')
hold on
title('Net Thrust vs Altitude vs Mach Number')
xlabel('Altitude ')
ylabel('Net Thrust / N ')

figure(2)
plot(net_thrust_array(i,:), sfc_array(i,:) ,'-o')
hold on
title('Specific Fuel Consumption vs Net Thrust')
xlabel('Net Thrust / N ')
ylabel('Specific Fuel Comsumption (SFC) / kg/s/N ')

figure(3)
plot(height_array(i,:) , spec_net_thrust_array(i,:) ,'-o')
hold on
title('Specific Net Thrust vs Altitude vs Mach Number')
xlabel('Altitude ')
ylabel('Spec Net Thrust / N.s/kg ')

figure(4)
plot(height_array(i,:) , sfc_array(i,:) ,'-o')
hold on
title('Specific Fuel Comsumption vs Altitude vs Mach Number')
xlabel('Altitude ')
ylabel('Specific Fuel Comsumption (SFC) / kg/s/N ')

figure(5)
plot(height_array(i,:) , air_mass_array(i,:) ,'-o')
hold on
title('Inlet Air Flow vs Altitude vs Mach Number')
xlabel('Altitude ')
ylabel('Inlet Air Mass Flowrate / kg/s ')

i = i+1;


end

sfc_trans = sfc_array.';
net_thrust_trans = net_thrust_array.';

for x = 1:num
    for y = 1:j
    figure(2)
plot(net_thrust_trans(x,:), sfc_trans(x,:) ,'-o')
hold on    
    end
end

end


%% Functions attached

function [T2 , P2, W] = Comp_Func (  T1 , P1 ,PR , Eff , gamma , Cpa )
    % Calculates the massflow , temperature and pressure of a compressor with a
    % given massflow , pressure , pressure ratio , efficiency and using the gamma
    % and CP of air.
    % Calculate Pressure
    P2=PR*P1;
    % Calculate Temperature using Isentropic efficiency
    T2 = T1*(1/Eff)*((P2/P1)^(( gamma -1)/ gamma ) -1) + T1;
    % Calculate work done by the compressor
    W=Cpa *(T2 -T1 );
end
    function [ T2 ,P2 , AF ] = burner_func ( T1 , P1 , Cpg, PR_burn, BET, Qhv )

T2=BET;
P2=PR_burn*P1;
%AF
delta_T = T2-T1;
AF = (Qhv/(Cpg * delta_T))-1;
    end
    
function [ m2 , T2 , P2 , Velo_2 ] = Nozzle_Func( T1, P1, Eff_noz, gamma_gas, A_noz,R )
% Calculates the massflow , temperature and pressure of an exhaust with a given
%  pressure , pressure ratio , efficiency and using the gamma of gas,
%  Area,entropy at inlet

%Finding Pc
exponent = (gamma_gas)/(gamma_gas-1);
ratio = (1-(1/Eff_noz)*((gamma_gas-1)/(gamma_gas+1)))^exponent ;
Pc = P1*ratio;
%Calculate T2
T2 = T1*(2*(1/(gamma_gas + 1)));

%Calculating exit pressure
P2= P1*(ratio);

%calcuate density
row = Pc/(R*T2);
%calcuate velocity
Velo_2 = (gamma_gas*R*T2)^0.5;
%Calculating mass flow rate
m2 = A_noz * row * Velo_2;

end

function [ T2 , P2 ] = Turb_Func( T1 , P1 , Eff_isen , Eff_mechn , W , gamma_gas , Cpg )
% Calculates the massflow , temperature and pressure of a turbine with a given
% massflow , pressure , Mechanical and isentropic efficiencies , the work done ,
% the cp of gass and using the gamma of gas
% Calculate temperature difference due to work done
Delta_T = W/(Cpg*Eff_mechn);
%Calculating T2
T2 = T1 - Delta_T;
%Calculate P2 with isentropic effeciency
P2 = (1/Eff_isen)*P1*((T2/T1)^(gamma_gas/(gamma_gas-1)) - 1 ) + P1;
end

function [thrust_net , thrust_s] = thrust_func(Velo_1,Velo_2, m2,Po,p2,A_noz,m_inlet)

%Calculating Thrust

%gross thrust

pressure_thrust = (A_noz)*(p2-Po);

thrust_g = m2*Velo_2;

%Ram Drag

ram_drag = m_inlet*Velo_1;

%Net Thrust

thrust_net = thrust_g + pressure_thrust - ram_drag ;

%Spec Net Thrust

thrust_s = thrust_net/m_inlet;



end