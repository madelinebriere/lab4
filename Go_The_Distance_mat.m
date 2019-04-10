clear

%% Time Vectors

fs=8000; % sampling frequency
whole = (1-1/fs)*2.4; %approximate length of measure
t1 = [0:1/fs:whole+2/fs]; %whole note - adjusted to have an easy-to-work with vector length
t2 = [0:1/fs:whole/2+1/fs]; %half note
t2dot = [0:1/fs:(whole)*(3/4)+1/fs]; %dotted half note
t3 = [0:1/fs:whole*(3/4)+1/fs]; 
t4 = [0:1/fs:whole/4]; %quarter notes
t4dot = [0:1/fs:whole*(3/8)]; %dotted quarter note
ttrip = [0:1/fs:whole/12]; %triplet
t8 = [0:1/fs:whole/8]; %eighth note
t8dot = [0:1/fs:whole*(3/16)]; %dotted eighth note
t16 = [0:1/fs:whole/16]; %sixteenth note
t16dot = [0:1/fs:whole*(5/16)]; %dotted sixteenth note
rest8 = zeros(1,length(t8)); %eighth rest
rest16 = zeros(1,length(t16)); %sixteenth rest
bet = [0:1/fs:.001]; %length of pause b/w identical notes
bw = zeros(1,length(bet)); %zero vector for pause

%% Note Frequencies

A = 220;
As = 220*2^(2/12);
Ah = A*2;
Al = A/2;
B = 220*2^(2/12);
Bl = B/2;
C = 220*2^(3/12);
Cs = 220*2^(4/12);
Csl = Cs/2;
D = 220*2^(5/12);
Dl = D/2;
E = 220*2^(7/12);
El = E/2;
F = 220*2^(8/12);
Fs = 220*2^(9/12);
Fsl = Fs/2;
G = 220*2^(10/12);
Gs = 220*2^(11/12);
Gl = G/2;
Af = Gs/2;


%% Song

% Melody function using harmonics, parameters of frequency and time
m = @(f,time) (0.1155*cos(2*pi*f*time+2.1299)+.3417*cos(2*pi*2*f*time-1.6727) ...
    + .1789*cos(2*pi*3*f*time+2.5454)+.1232*cos(2*pi*4*f*time-0.6607) ...
    +.0678*cos(2*pi*5*f*time+2.0390)+.0473*cos(2*pi*f*6*time-2.1597) ...
    +.0260*cos(2*pi*f*7*time+1.0467)+.0045*cos(2*pi*f*8*time-1.8581) ...
    +.002*cos(2*pi*f*9*time+2.3925)).*(exp(-7.*time).*(6.*time.^2+5.*time));
%https://web.eecs.umich.edu/~fessler/course/100/misc/course-notes-ay-jf.pdf

% Harmony function using instrument synthesis
h = @(f,t) 0.025*exp(-0.6*t).*sin(2*pi*f*t+3.5*exp(-0.6*t).*sin(2*pi*f*t*2));
% Harmony function with faster fade used for end
hend = @(f,t) 0.025*exp(-2*t).*sin(2*pi*f*t+3.5*exp(-2*t).*sin(2*pi*f*t*2));

%% Melody, number after m corresponds to measure number
% Two bw per measure for consistent vector lengths
m1 = [rest8 rest8 rest8 rest8 rest8 rest8 rest8 bw bw m(Fsl,t16) m(A,t16)];
m2 = [m(B,t4) m(A,t4) m(Fsl,t4) rest8 m(Fsl,t16) m(A,t16) bw bw];
m3 = [m(B,t4) m(A,t4) m(Fsl,t4) rest8 m(Fsl,t16) m(A,t16) bw bw];
m4 = [m(B,t4) m(Cs,t4) m(D,t8) m(A,t8dot) rest16 m(Dl,t16) m(Fsl,t16) bw bw];
m5 = [bw m(Fsl,t4) bw m(Fsl,ttrip) m(Gl,ttrip) m(Fsl,ttrip) m(El,t4) rest8 m(Fsl,t16) m(A,t16)];
m6 = [ m(B,t4) m(A,t16) bw m(A,t8) bw m(A,t16) m(Fsl,t4dot) m(D,t16) m(Cs,t16)];
m7 = [m(B,t4) m(A,t8) bw m(A,t16) m(Fsl,t16) m(Dl,t4dot) m(Fsl,t16) m(A,t16) bw];
m8 = [m(B,t4) m(Cs,t4) m(D,t16) m(B,t8dot) m(D,t16) m(Fs,t8dot) bw bw];
m9 = [m(Fs,t8) m(G,t8) m(Fs,t8) m(G,t16) m(E,t16) m(E,t4) rest8 m(Fs,t16) m(G,t16) bw bw];
m10 = [0.8*m(Ah,t4) m(D,t8) m(E,t2) rest8 bw bw];
m11 = [m(Fs,t8) m(E,t16) m(Fs,t8dot) m(G,t16) m(Fs,t8), m(E,t8dot) rest8 m(Fs,t16) m(G,t16) bw bw];
%High A note is too high of a frequency, so reduced the amplitude
m12 = [0.8*m(Ah,t4) m(D,t8) m(E,t2) rest8 bw bw];
m13 = [m(Fs,t16) m(E,t8) m(Fs,t8) m(G,t16) m(Fs,t8) m(E,t4) rest8 m(Fs,t16) m(G,t16) bw bw];
m14 = [0.8*m(Ah,t8) m(D,t8) m(Cs,t8) m(B,t2) m(D,t16) m(E,t16) bw bw];
m15 = [m(Fs,t4) m(A,t8) bw m(A,t2) rest8 bw];
m16 = [rest8 rest8 rest8 rest8 rest8 rest8 rest8 bw bw m(Fsl,t16) m(A,t16)];
m17 = [m(B,t4) m(A,t4) m(Fsl,t8) m(A,t8) m(D,t8) m(Fs,t8) bw bw];
m18 = [m(Fs,t4dot) m(G,t8) m(E,t8) m(D,t16) m(B,t8dot) m(0,t16) bw bw m(D,t16)];
m19 = [m(D,t1) bw bw];

% Melody + reverb
mo = [m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19 rest16];
mr = 0.15.*([rest16 m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19]);

% First harmony. Variable name is h1 + measure number (corresponds to
% melody)
h11 = [h(0,t1) bw bw];
h12 = [bw bw h(Gl,t4) h(El,t4) h(Dl,t2)];
h13 = [h(Gl,t4) h(El,t4) h(Dl,t2) bw bw];
h14 = [h(Gl,t4) h(A,t4) h(Fsl,t4) rest8 h(Dl,t16) h(Fsl,t16) bw bw];
h15 = [bw h(Dl,t2) bw h(El,t4dot) h(Fsl,t16) h(A,t16)];
h16 = [h(Gl,t4) h(El,t4) bw bw h(Dl,t2)];
h17 = [h(Gl,t4) h(El,t4) bw rest8 rest8 rest8 h(Fsl,t16) h(A,t16) bw];
h18 = [h(Gl,t4) h(As,t4) h(Fsl,t4) h(As,t4) bw bw];
h19 = [h(D,t3) rest8 h(Fs,t16) h(G,t16) bw bw];
h110 = [0.8*h(E,t4) h(A,t4) h(B,t4) h(Cs,t8) h(D,t8) bw bw];
h111 = [h(D,t2) h(Cs,t2) bw bw];
h112 = h110;
h113 = h111;
h114 = [0.8*h(D,t4) h(Cs,t8) h(B,t8) h(Gl,t2) bw bw];
h115 = [h(D,t4) h(A,t8) bw h(A,t2) rest8 bw];
h116 = [bw bw hend(Dl,t2dot) h(Dl,t8) h(Fsl,t16) h(A,t16)];
h117 = [h(Gl,t4) h(El,t4) h(Dl,t8) h(Dl,t8) h(B,t8) h(D,t8) bw bw];
h118 = [h(B,t2) h(B,t2) bw bw];
h119 = [h(A,t4) h(A,t4) hend(Gl,t2) bw bw rest16];

% Second harmony. Same syntax as first.
h21 = [h(0,t1) bw bw];
h22 = [bw bw h(Dl,t4) h(Csl,t4) h(Al,t2)];
h23 = [h(Gl,t4) h(El,t4) h(Dl,t2) bw bw];
h24 = [h(Dl,t4) h(El,t4) h(Dl,t4) rest8 h(Dl,t16) h(Fsl,t16) bw bw];
h25 = [bw h(Bl,t2) bw rest8 rest8 h(Csl,t4)];
h26 = [h(Dl,t4) h(Csl,t4) bw bw h(Al,t2)];
h27 = [h(Dl,t4) h(Csl,t4) bw h(Dl,t4) h(Bl,t4) bw];
h28 = [h(Dl,t4) h(El,t4) h(Dl,t4) h(Fsl,t4) bw bw];
h29 = [h(B,t3) h(Cs,t4) bw bw];
h210 = [0.8*h(A,t4) h(D,t8) h(E,t8) h(El,t2) bw bw];
h211 = [h(A,t2) h(A,t2) bw bw];
h212 = h210;
h213 = h211;
h214 = [0.8*h(A,t4) h(Dl,t4) h(Dl,t2) bw bw];
h215 = [h(A,t4) h(Dl,t2dot) bw bw];
h216 = [bw bw hend(Gl,t2dot) h(Gl,t8) rest8];
h217 = [h(Dl,t4) h(Csl,t4) h(Al,t8) h(Csl,t8) h(Fsl,t8) h(A,t8) bw bw];
h218 = [h(Gl,t2) h(Gl,t2) bw bw];
h219 = [h(Dl,t4) h(Csl,t4) hend(Bl,t2) bw bw rest16];

% Concatenation of melody and each of the harmonies
allm = [m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19];
allh1 = [h11 h12 h13 h14 h15 h16 h17 h18 h19 h110 h111 h112 h113 h114 h115 ...
    h116 h117 h118 h119];
allh2 = [h21 h22 h23 h24 h25 h26 h27 h28 h29 h210 h211 h212 h213 h214 h215 ...
    h216 h217 h218 h219];

%soundsc(allh1 + allh2 + mo + mr); %for playing all together

%% .Wav File
% 
% %final vector resulting from adding all together
% final = mo + mr + allh1 + allh2;
% filename = 'Horowitz_GoTheDistance.wav';
% y = final./max(final);
% audiowrite(filename,y,fs);
% 
% %% Plots
% % All of the melody
% clf
% figure(1)
% plot(t4,m(15,t4));
% xlabel('time')
% % 
% % % Just harmonics of the melody
% mfun1 = @(f,time) (0.1155*cos(2*pi*f*time+2.1299)+.3417*cos(2*pi*2*f*time-1.6727) ...
%     + .1789*cos(2*pi*3*f*time+2.5454)+.1232*cos(2*pi*4*f*time-0.6607) ...
%     +.0678*cos(2*pi*5*f*time+2.0390)+.0473*cos(2*pi*f*6*time-2.1597) ...
%     +.0260*cos(2*pi*f*7*time+1.0467)+.0045*cos(2*pi*f*8*time-1.8581) ...
%     +.002*cos(2*pi*f*9*time+2.3925));
% figure(2)
% plot(t4,mfun1(15,t4));
% xlabel('time')
% % 
% % % Just the window of the melody
% mfun2 = @(time) (exp(-7.*time).*(6.*time.^2+5.*time));
% figure(3)
% plot(t4,mfun2(t4));
% xlabel('time')
% 
% % % All of the harmony
% figure(4)
% plot(t4,h(15,t4))
% xlabel('time')
% % 
% % % Just synthesis of the harmony
% hfun1 = @(f,t) sin(2*pi*f*t+sin(2*pi*f*t*2)); 
% figure(5)
% plot(t4,hfun1(15,t4)) 
% xlabel('time')
% % Just the window of the harmony 
% hfun2 = @(t) 0.025*exp(-0.6*t); 
% figure(6)
% plot(t4,hfun2(t4));
% xlabel('time')