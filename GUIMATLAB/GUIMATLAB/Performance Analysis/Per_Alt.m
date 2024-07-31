function [net_thrust_array , sfc_array] = Per_Alt(alt_1 , alt_2 , mach, comp_ratio)
num = (alt_2-alt_1)/1000 + 1;
Temp = zeros(5,1);
Press = zeros(5,1);
net_thrust_array = zeros(1,num);
spec_net_thrust_array = zeros(1,num);
sfc_array = zeros(1,num);
Alt_array = zeros(1,num);
air_mass_array = zeros(1,num);


%%[R,gamma_air,gamma_gas,cp_air,cp_gas,hf]= standard_vals;
%%[ m2 , T2 , p2 , W] = Comp_Func ( Mo , To , Po ,pi_c , Eff_c , gamma_gas , cp_gas );

cp_gas = 1000;  %J/Kg.Kelvin
Qhv = 43260000; %J/kg
gamma_gas = 1.4;
R = 287;
Eff_c = 1;
BET=1050;
PR_burn = 0.96;
Eff_Noz = 1;
A_noz = 0.3;
Eff_isen_turb = .9;
Eff_mechn_turb = 1;

i = 1;
j = 1;

%[Po,To,Mo,cp_gas,tau_r,tau_d,tau_c,tau_b,tau_t,tau_n,pi_r,pi_d,pi_c,pi_b,pi_t,pi_n,Ttb_lim,Eff_c,Eff_b,Eff_t,Eff_n] = inputparameters ;

%[Po,To,Mo,mo,cp_gas, gamma_gas] = [100e+3 , 223 , 2 , 32 , 1000 , 1.4];





j = 1;

for altitude = alt_1:1000:alt_2

    h = altitude;
    Tamb = (15.04-(0.00649*h));
    Pamb = (101.29* ((Tamb+273.1)/288.08)^5.256);
Po = Pamb*1000; %Pascals
To = Tamb + 273.16; %Kelvin
Velo_sound = sqrt(gamma_gas*R*(Tamb + 273));
Temp(1,1) = To;
Press(1,1) = Po;

Mo = mach;
PR_c = comp_ratio;

Alt_array(i,j) = altitude;
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
m_inlet = m5-m_fuel;
air_mass_array(i,j) = m_inlet;
specific_fuel = m_fuel/net_thrust;
sfc_array(i,j) = specific_fuel;



j = j + 1;

end



%PLOTS

figure(1)
plot(Alt_array(i,:), net_thrust_array(i,:) ,'-o')
hold on
title('Net Thrust vs Altitude')
xlabel('Altitude / m ')
ylabel('Net Thrust / N ')

figure(2)
plot(net_thrust_array(i,:), sfc_array(i,:) ,'-o')
hold on
title('Specific Fuel Consumption vs Net Thrust')
xlabel('Net Thrust / N ')
ylabel('Specific Fuel Comsumption (SFC) / kg/s/N ')

figure(3)
plot(Alt_array(i,:) , spec_net_thrust_array(i,:) ,'-o')
hold on
title('Specific Net Thrust vs Altitude')
xlabel('Altitude / m  ')
ylabel('Spec Net Thrust / N.s/kg ')

figure(4)
plot(Alt_array(i,:) , sfc_array(i,:) ,'-o')
hold on
title('Specific Fuel Comsumption vs Altitude')
xlabel('Altitude / m  ')
ylabel('Specific Fuel Comsumption (SFC) / kg/s/N ')

figure(5)
plot(Alt_array(i,:) , air_mass_array(i,:) ,'-o')
hold on
title('Inlet Air Flow vs Altitude')
xlabel('Altitude / m ')
ylabel('Inlet Air Mass Flowrate / kg/s ')




end


%% Fucntions Attached

    
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