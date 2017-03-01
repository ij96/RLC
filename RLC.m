function RLC(Vin_str,method)
    % Plot RLC series circuit response by using numerical method
    %
    % The RLC series circuit can be described by the equation below:
    %   L*q'' + R*q' + (1/C)*q = Vin
    %   where:
    %       q(t): charge
    %       Vin(t): input voltage
    %       R,L,C: values of R (in Ohms), L (in Henrys) and C (in Farads)
    %
    %   The equation above can be written as two coupled 1st order ODEs:
    %   q' = i
    %   i' = (1/L)*(Vin - R*i - (1/C)*q)
    %   where:
    %       i(t): current = dq/dt = q'
    %
    % Parameters:
    %   Vin_str: Vin(t) as a string
    %   method: name of the RK method to be used (case-insensitive)
    %
    % Add-ons required:
    %   Symbolic Math Toolbox

    % Command line argument
    syms t;                     % create a symbol, needed to parse Vin_str
    s = eval(Vin_str);          % converts string to epression/constant
    s_type = char(class(s));    % get type for if condition
    Vin_is_double = strcmp(s_type,'double');

    % Convert Vin_str into Vin, either as a const or function of t
    if(Vin_is_double)
       Vin = s;                 % Vin as a const
    else
       Vin = matlabFunction(s); % Vin as a function of t
    end

    % Circuit parameters
    R = 260;                    % Ohm
    L = 0.5;                    % Henry
    C = 0.0000035;              % Farad

    % Circuit initial conditions
    qi = 0.000000005;           % Coulomb
    ii = 0;                     % Ampere

    % Set time range and interval size
    ti = 0;                     % initial time
    tf = 0.05;                  % final time
    h = 0.0001;                 % time interval size

    % Coupled 1st order ODEs
    dq = @(t,q,i) i;
    if(Vin_is_double)
       di = @(t,q,i) (1/L)*(Vin - R*i - (1/C)*q);
    else
       di = @(t,q,i) (1/L)*(Vin(t) - R*i - (1/C)*q);
    end

    % Applying numerical method
    [t,q,i] = RK(method,dq,di,qi,ii,ti,tf,h);

    % Plot graph
    V_C = (1/C)*q;
    V_R = R*i;
    figure;
    if(Vin_is_double)
        plot(t, ones(size(t))*Vin, 'k:');
        V_L = Vin - V_R - V_C;
    else
        plot(t, Vin(t),'k:');
        V_L = Vin(t) - V_R - V_C;
    end
    hold on;
    plot(t, V_C, t, V_R, t, V_L);

    title('Voltage across R, L and C');
    xlabel('Time/s');
    ylabel('Voltage/V');
    legend('Vin(t)','V_C = (1/C)*q','V_R = R*q''','V_L = L*q''''');