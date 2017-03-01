function s = Square(t,freq)
    % Custom square wave function
    % duty cycle: 50%, amplitude: 0 to 1
    % Use this instead of 'square', since Syms does not work with the 'square' func in Signal Processing Toolbox
    % 
    % t: time
    % freq: in Hertz
    s = 2*floor(freq*t)-floor(2*freq*t)+1;