function RL(Vin_str,method)
    % Plot RL series circuit response by using numerical method
    %
    % The RL series circuit can be described by the equation below:
    %   L*i' + R*i = Vin
    %   where:
    %       i(t): current
    %       Vin(t): input voltage
    %       R,L: values of R (in Ohms) and L (in Henrys)
    %
    % Parameters:
    %   Vin_str: Vin(t) as a string
    %   method: name of the RK method to be used (case-insensitive)
    %
    % Add-ons required:
    %   Symbolic Math Toolbox
    %
    % Note:
    %   This function uses RK.m which is written for solving 2nd order
    %   ODE. Two dummy variables (dq and qi) are introduced as the
    %   parameters required by RK.m

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

    % Circuit initial conditions
    qi = 0;                     % dummy initial condition
    ii = 0;                     % Ampere

    % Set time range and interval size
    ti = 0;                     % initial time
    tf = 0.05;                  % final time
    h = 0.0001;                 % time interval size

    % Coupled 1st order ODEs
    dq = @(t,q,i) 1;            % dummy 1st order ODE
    if(Vin_is_double)
       di = @(t,q,i) (1/L)*(Vin - R*i);
    else
       di = @(t,q,i) (1/L)*(Vin(t) - R*i);
    end

    % Applying numerical method
    [t,q,i] = RK(method,dq,di,qi,ii,ti,tf,h);

    % Plot graph
    V_R = R*i;
    figure;
    if(Vin_is_double)
        plot(t, ones(size(t))*Vin, 'k:');
        V_L = Vin - V_R;
    else
        plot(t, Vin(t),'k:');
        V_L = Vin(t) - V_R;
    end
    hold on;
    plot(t, V_R, t, V_L);

    title('Voltage across R and L');
    xlabel('Time/s');
    ylabel('Voltage/V');
    legend('Vin(t)','V_R = R*i','V_L = L*i''');