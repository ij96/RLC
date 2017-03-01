function RLC(Vin_str)
% Plot RLC series circuit response by using numerical method (Runge-Kutta 3/8)
%
% RLC system function:
%   L*q'' + R*q' + (1/C)*q = Vin
%   where:
%       q(t): charge
%       Vin(t): input voltage
%
%   Write the equation above as coupled 1st order ODEs:
%   q' = i
%   i' = (1/L)*(Vin - R*i - (1/C)*q)
%   where:
%       i(t): current = dq/dt = q'
%
% Add-ons required:
%   Symbolic Math Toolbox

% command line argument
syms t;                     % create a symbol, needed to parse Vin_str
s = eval(Vin_str);          % converts string to epression/constant
s_type = char(class(s));    % get type for if condition

% convert Vin_str into Vin, either as a const or function of t
if(strcmp(s_type, 'double'))
   Vin = s;                 % Vin as a const
else
   Vin = matlabFunction(s); % Vin as a function of t
end

% circuit parameters
R = 260;                    % Ohm
L = 0.5;                    % Henry
C = 0.0000035;              % Farad

% circuit initial conditions
qi = 0.000000005;           % Coulomb
ii = 0;                     % Ampere

% set time range and interval size
ti = 0;                     % initial time
tf = 0.05;                  % final time
h = 0.0001;                 % interval size

% coupled 1st order ODEs
dq = @(t,q,i) i;
if(strcmp(s_type, 'double'))
   di = @(t,q,i) (1/L)*(Vin - R*i - (1/C)*q);
else
   di = @(t,q,i) (1/L)*(Vin(t) - R*i - (1/C)*q);
end

% applying numerical method
N = round((tf-ti)/h);       % total number of steps
t = zeros(1,N);             % [1 x N] matrix of 0's
q = zeros(1,N);
i = zeros(1,N);

% set initial conditions
t(1) = ti;
q(1) = qi;
i(1) = ii;

% iterate through N step and get the [i+1]th value each iteration.
for j = 1:N-1
    [q(j+1),i(j+1)] = RK4_38(dq,di,t(j),q(j),i(j),h);
    t(j+1) = t(j) + h;
end

% plot graph
V_C = (1/C)*q;
V_R = R*i;
figure;
if(strcmp(s_type, 'double'))
    plot(t, ones(size(t))*Vin,'k:');
else
    plot(t, Vin(t),'k:');
end
hold on;
if(strcmp(s_type, 'double'))
   V_L = Vin - V_R - V_C;
else
   V_L = Vin(t) - V_R - V_C;
end
plot(t, V_C, t, V_R, t, V_L);

title('Voltage across R, L and C');
xlabel('Time/s');
ylabel('Voltage/V');
legend('Vin(t)','V_C = (1/C)*q','V_R = R*q''','V_L = L*q''''');