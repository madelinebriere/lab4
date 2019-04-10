Fs = 44100;
T = 1/Fs;
N = 400;
f = Fs/N; % 44100 / 400 ~ A

fr = zeros(1, N/2);
fl = zeros(1, N/2);

% Sample test waveform.
% start = sin(linspace(0, pi, N/2))/2;
% fr = start;
% Piecewise waveform.
pluck_point = N/10;
sample_point = N/4;
start = [linspace(0, 1, pluck_point) linspace(1, 0, N/2-pluck_point)];
fr = start;
fl = start;

pt = .01;
l = 2;
record_length = l * Fs;
recorded = zeros(1, record_length);
for i = 1:record_length
    % pause(pt)
    % plot(x, fr + fl)
    % xlim([0 N/2])
    % ylim([-1 1])
    x = 1:N/2;
    recorded(1, i) = fr(1, sample_point) + fl(1, sample_point);
    fl_temp = fl;
    fr_temp = fr;
    fl(1, N/2) = -fr_temp(1, N/2);
    fr(1, 1) = -fl_temp(1, 1);
    fr(1, 2:N/2-1) = fr_temp(1, 1:N/2-2);
    fl(1, 1:N/2-1) = fl_temp(1, 2:N/2);
    % Filter.
    a = .02;
    y = .99;
    fr(1, N/2) = y.*(-a.*fl_temp(1, N/2) + (1-a).*fr_temp(1, N/2-1));
end 
soundsc(recorded, Fs);
