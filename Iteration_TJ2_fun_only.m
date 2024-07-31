function Iteration_TJ2_fun_only

clear,clc,close all;

Temp = zeros(6,1);
Press = zeros(6,1);
s_vol = zeros(6,1);
s_entro = zeros(6,1);

%[Po,To,Mo,cp_gas,tau_r,tau_d,tau_c,tau_b,tau_t,tau_n,pi_r,pi_d,pi_c,pi_b,pi_t,pi_n,Ttb_lim,Eff_c,Eff_b,Eff_t,Eff_n] = inputparameters ;

%[Po,To,Mo,mo,cp_gas, gamma_gas] = [100e+3 , 223 , 2 , 32 , 1000 , 1.4];

cp_gas = 1000;  %J/Kg.Kelvin
Qhv = 43260000; %J/kg
gamma_gas = 1.4;
R = 287;
h = 10000;
Mo = .6;
Tamb = (15.04-(0.00649*h));
Pamb = (101.29* ((Tamb+273.1)/288.08)^5.256);
Velo_sound = (gamma_gas*R*(Tamb + 273))^0.5;
Velo_1 = Mo*Velo_sound;
Dynamic_Temp = (Velo_1^2)/(2*cp_gas);
Tstp = 273.16;
Pstp = 101325;
Ndp = 50700;

Pt = Pamb*1000*(1 + Dynamic_Temp/(Tamb + 273))^(gamma_gas/(gamma_gas-1)); %Pascals
Tt = Tamb + 273.16 + Dynamic_Temp; %Kelvin
Po = Pt; %Pascals
To = Tt; %Kelvin



% so  = 1000*(air_table(Tt)); % This needs to be automated
so = 1340;

Temp(1,1) = To;
Press(1,1) = Po;
s_entro(1,1) = so;
vo=(R*To)/Po;
s_vol(1,1) = vo;

T1 = To;
p1 = Po;
s1 = so;
v1 = vo;

%Compressor Function Call
Temp(2,1) = T1;
Press(2,1) = p1;
s_entro(2,1) = s1;
s_vol(2,1) = v1;

Eff_c = .85;
PR_c = 8; 
[ T2 , p2 ,v2 , s2 , W_c ] = Comp_Func (T1 , p1 , PR_c , Eff_c , gamma_gas , cp_gas , R , s1);


Temp(3,1) = T2;
Press(3,1) = p2;
s_entro(3,1) = s2;
s_vol(3,1) = v2;

%Burner Function Call

BET=1050;
PR_burn = 0.96;
[ T3 ,p3 , v3 , s3,AF] = burner_func (  T2 , p2 ,  cp_gas , R , s2,PR_burn,BET,Qhv );



Temp(4,1) = T3;
Press(4,1) = p3;
s_entro(4,1) = s3;
s_vol(4,1) = v3;

%cp_gas = 1004;

%Turbine Function Call

Eff_isen_turb = .9;
Eff_mechn_turb = 1;
W_t = W_c;

[ T4 , p4 ,v4 , s4 ] = Turb_Func( T3 , p3 , Eff_isen_turb , Eff_mechn_turb , W_t , gamma_gas , cp_gas , s3 ,R);



Temp(5,1) = T4;
Press(5,1) = p4;
s_entro(5,1) = s4;
s_vol(5,1) = v4;

%Nozzle Function Call

Eff_Noz = 1;
%PR_noz = .05;
A_noz = 0.06;

[ m5 , T5 , p5 , Velo_5 , v5 , s5 ] = Nozzle_Func( T4 , p4 , Eff_Noz , gamma_gas , cp_gas , s4, R , A_noz);

dyn_temp = ((Velo_5)^2)/(2 * cp_gas);
dyn_press = ((Velo_5)^2)/(v5*2);
Temp(6,1) = T5 + dyn_temp;
Press(6,1) = p5  + dyn_press;
s_entro(6,1) = s5;
s_vol(6,1) = v5;


%Thrust Function Call
m_fuel = m5 /(1+AF);
m_inlet = m5-m_fuel;

[net_thrust , net_thrust_spec] = thrust_func(Velo_1,Velo_5, m5,Po,p5,A_noz,m_inlet);
thrust1 = net_thrust;

%SFC Calculation

specific_fuel = m_fuel/net_thrust;
sfc1 = specific_fuel;

% Other cals
pie_dp = PR_c;
eta_t = 0.9;
gamma = 1.4;%%HUZAIFAH SEEEE
a=1.5;
b=5;
k=0.03;
c=3;
d=4;
C_eff=15;
D=1;
eeta_o = .85;
theeta1 = T1/Tstp;
theeta3=T3/Tstp;
theeta5=T5/Tstp;
mcor_1 =((m_inlet * ((theeta1)^0.5))/(p1/Pstp));
mcor_3 =((m_inlet * ((theeta3)^0.5))/(p3/Pstp));
mcor_5 =((m_inlet * ((theeta5)^0.5))/(p5/Pstp));
pic = p2/p1;
pib = p3/p2;
pit=p4/p3;
Ncor_c = Ndp/((theeta1)^0.5);
Ncor_t = Ndp/((theeta3)^0.5);
A6 = mcor_5/mcor_1;
B6 =(((T3/T1) * (T4/T3)) * abs(1/(pic*pib*pit)))^0.5;
C6 = ((T4/T1)*(1/(p4/p1)))^0.5;
mcor_4NR =mcor_3*Ncor_t;
mcor_4NR2 =mcor_1 *(1/(pic*pib))*Ncor_c;

rpm_array = [0.76 0.82 0.875 0.922 0.964 1 1.033 1.063 ];

rows = size(rpm_array,1);
cols = size(rpm_array,2);
total_array = zeros(cols , 4);
hit = 0;

for x = 5 : 5
    
  pi_max = (1+(pie_dp-1)*(((rpm_array(x))^(a*b))+2*rpm_array(x)*k*log(1-((0-(rpm_array(x)^b))/k))));
    
Convergence1 = 0;

guess=(pi_max)-.1;

% guess=pi_max-.1;
Ncorrcg1 = rpm_array(x)*Ncor_c;
Ncg1 = Ncorrcg1 * (theeta1)^0.5 ;
aaa = 0;
picg1 = 1;
p_hat = 0.5;

Convergence_array = zeros(50,1);
m_hat_array = zeros(50,1);
eeta_array = zeros(50,1);
p_hat_array = zeros(50,1);
count = 0;

while aaa == 0 && picg1 >= 0
    
count = count + 1;    
picg1=guess;
delta_a=0;

syms Xmass
p_hat = (picg1 -1)/(pie_dp-1);

m_hat = solve((1+(pie_dp-1)*(((rpm_array(x))^(a*b))+2*rpm_array(x)*k*log(1-((Xmass-(rpm_array(x)^b))/k))))-picg1 , Xmass);

mo_hat = m_hat;

eeta_req = eeta_o * (1 - (C_eff* (abs(( (p_hat/(m_hat^(a+delta_a-1))) -m_hat)))^c ) - D*((abs((m_hat/mo_hat)-1))^d)     );

etacg1= eeta_req;
taucg1 = (((picg1^((gamma-1)/(gamma)))-1)/etacg1)+1;
mdotcorr2g1 = m_hat*mcor_1;
mdotcorr4normNg1 = mdotcorr2g1*(1/(picg1*pib))*Ncorrcg1;
% mdotcorr4normNgpercent1 = (mdotcorr4normNg1/mcor_4NR)*100;
Ncorrtg1 = mdotcorr4normNg1/(mcor_3);
Tt4Tt2g1 = (Ncorrcg1/Ncorrtg1)^2;
tautg1 = 1 - (1/(Tt4Tt2g1))*(taucg1-1);
pitg1 = (1-(1-tautg1)/(eta_t))^(gamma/(gamma-1));

%matching criteria
A1 = mcor_5/mdotcorr2g1;
B1 = ((Tt4Tt2g1*tautg1)^0.5)*abs(1/(picg1*pib*pitg1));
AA = round(A1,7);
BB = round(B1,7);

% if x == 1
%     toll = 0.00005;
%     
% elseif x>=2 && x<=3
%         toll = 0.0005;
% else
%         toll = 0.005;
% end
% 

Convergence1 = (abs(A1-B1));
guess = guess - 0.005;
ratio = picg1/pie_dp;
% 

Y = round(Convergence1,5);
z = round(m_hat,5);
W = round(eeta_req,5);


% if Y > 100
% 
%     hit = 1;
%     
% end


% if x == 1
%     conver_criteria = 0.1;
% else
%     conver_criteria = 0.05;
% end
% 
% 
% if hit==1
%     
%     if Y <= conver_criteria
%         
%     aaa = 1;
%     total_array(x,1) = Ncg1;
%     total_array(x,2) = ratio;
%     total_array(x,3) = m_hat;
%     total_array(x,4) = etacg1;
%     fprintf('The value of m_hat is %.3f\n', m_hat);
%     fprintf('The value of etacg1 is %.3f\n', etacg1);
%         
%     end
%     
% end




% Convergence_array(count,1) = Y;
% m_hat_array(count,1) = z;
% eeta_array(count,1) = W;
% p_hat_array(count,1) = p_hat;
% 
% [AA,BB]=meshgrid(min(m_hat_array):.01:max(m_hat_array),min(p_hat_array):.01:max(p_hat_array));

% xg = linspace(min(m_hat_array), max(m_hat_array), 5);
% yg = linspace(min(p_hat_array), max(p_hat_array), 5);
% [AA,BB]=meshgrid(xg,yg);

% CC = griddata(m_hat_array, p_hat_array, eeta_array, AA, BB);
% 

figure(1)
hold ON
plot(z,ratio,'b*')



% if Y >= 0.05 && (AA>=BB)
%     guess = guess - 0.005;
% elseif Y >= 0.05 && (AA<=BB)
%     guess = guess + 0.005;
% else
%     aaa = 1;
%     total_array(x,1) = Ncg1;
%     total_array(x,2) = guess;
%     total_array(x,3) = mdotcorr2g1;
%     total_array(x,4) = etacg1;
%     fprintf('The value of m_hat is %.3f\n', m_hat);
% fprintf('The value of etacg1 is %.3f\n', etacg1);  
%     
% end

fprintf('The value of Convergence is %.3f\n', Convergence1);

figure(2)
hold ON
plot(Y,ratio,'r*')

end

% figure(1)
% hold ON
% plot(total_array(x,3),total_array(x,2),'r*')
% fprintf('The value of Convergence is %.3f\n', Convergence1);
% fprintf('The value of m_hat is %.3f\n', m_hat);
% fprintf('The value of etacg1 is %.3f\n', etacg1);  
end





%PLOTS
figure(1)
% figure(2)



% figure(5)
% contour(AA,BB,CC)

% figure(3)
% plot(s_entro,Temp ,'-o') 
% grid on;
% title('Temp vs Specific Entropy')
% xlabel('Specific Entropy / (kJ/kg.K)')
% ylabel('Temperature / K')
% zlabel('Pressure / Pa')
% 
% figure(4)
% plot(s_vol,Press ,'-o')
% grid on;
% title('Pressure vs Specific Volume')
% xlabel('Specific Volume / (m^3/kg)')
% ylabel('Pressure / Pa')

end


% Attached functions

function [thrust_net , thrust_s ] = thrust_func(Velo_1,Velo_2, m2,Po,p2,A_noz,m_inlet)

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

% function [T2,P2,s2,v2] =inlet_func(Velo_1,cp_gas,T1,P1,eff_inlet,s1,R,gamma)
% %delta_t= (Velo_1^2)/(2*cp_gas);
% delta_t = 0;
% T2 = T1+delta_t;
% pratio = (1 + eff_inlet*(delta_t/T1) )^(gamma/(gamma-1));
% P2 = pratio * P1;
% v2= (R*T2)/P2;
% delta_s = cp_gas*log((T2/T1)) - R*log(P2/P1);
% s2 = s1 + delta_s;
% end

    function [T2 , P2 , v2 , s2 , W] = Comp_Func (T1 , P1 ,PR , Eff , gamma , Cpa , R, s1)
    % Calculates the massflow , temperature and pressure of a compressor with a
    % given massflow , pressure , pressure ratio , efficiency and using the gamma
    % and CP of air.
    % Calculate Pressure
    P2=PR*P1;
    % Calculate Temperature using Isentropic efficiency
    T2 = T1*(1/Eff)*((P2/P1)^(( gamma -1)/ gamma ) -1) + T1;
    % Calculate work done by the compressor
    W = Cpa *(T2 -T1 );
    % Calculate v1 at Compressor Outlet
    v2 = (R*T2)/P2;
    % Calculate s1 at Compressor Outlet
    delta_s = Cpa*log((T2/T1)) - R*log(P2/P1);
    s2 = s1 + delta_s;

    end
    
function [ T2 ,P2 , v2, s2, AF ] = burner_func ( T1 , P1 , Cpg, R, s1, PR_burn, BET, Qhv )
T2=BET;
P2=PR_burn*P1;
% Calculate v2 at Burner Outlet
v2 = (R*T2)/P2;
%AF
delta_T = T2-T1;
AF = (Qhv/(Cpg * delta_T))-1;
% Calculate s2 at Burner Outlet
Delta_s = (Cpg)*log(T2/T1) - (R)*log(P2/P1);
s2 = Delta_s + s1;
end

function [ m2 , T2 , P2 , Velo_2 , v2 , s2 ] = Nozzle_Func(  T1 , P1 , Eff_noz , gamma_gas , Cpg ,s1 ,R, A_noz )
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
% Calculate v2 at nozzle Outlet
v2 = (R*T2)/P2;
% Calculate s2 at nozzle Outlet
delta_s=Cpg*log((T2/T1))- R*log(P2/P1);
s2 = s1 + delta_s;
end

function [ T2 , P2 , v2, s2 , PR ] = Turb_Func( T1 , P1 , Eff_isen , Eff_mechn , W , gamma_gas , Cpg, s1, R)
% Calculates the massflow , temperature and pressure of a turbine with a given
% massflow , pressure , Mechanical and isentropic efficiencies , the work done ,
% the cp of gass and using the gamma of gas
% Calculate temperature difference due to work done

Delta_T = W/(Cpg*Eff_mechn);
%Calculating T2
T2 = T1 - Delta_T;
%Calculate P2 with isentropic effeciency
P2 = (1/Eff_isen)*P1*((T2/T1)^(gamma_gas/(gamma_gas-1)) - 1 ) + P1;
%Calculating Pressure ratio
PR = P1/P2;
% Calculate v2 at Turbine Outlet
v2 = (R*T2)/P2;
% Calculate s2 at Compressor Outlet
delta_s = Cpg*log(T1/T2) - R*log(P1/P2);
s2 = s1 + delta_s;
end
    
    