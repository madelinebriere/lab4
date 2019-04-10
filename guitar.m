function output = guitar(f, l)
Fs = 44100;
N = round(Fs/f);
N = round(N/2)*2;

fr = zeros(1, N/2);
fl = zeros(1, N/2);

% Piecewise waveform.
pluck_point = round(N/10);
sample_point = round(N/4);
start = [linspace(0, 1, pluck_point) linspace(1, 0, N/2-pluck_point)];
fr = start;
fl = start;

record_length = l * Fs;
output = zeros(1, record_length);
for i = 1:record_length
    output(1, i) = fr(1, sample_point) + fl(1, sample_point);
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
end
