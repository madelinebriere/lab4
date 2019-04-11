fs=8000; % sampling frequency
whole = (1-1/fs)*2.4; %approximate length of measure
t1 = 2.4; %[0:1/fs:whole+2/fs]; %whole note - adjusted to have an easy-to-work with vector length
t2 = t1/2; %[0:1/fs:whole/2+1/fs]; %half note
t2dot = t1*3/4;%[0:1/fs:(whole)*(3/4)+1/fs]; %dotted half note
%t3 = [0:1/fs:whole*(3/4)+1/fs]; 
t4 = t1/4; %[0:1/fs:whole/4]; %quarter notes
t4dot = t1*3/8; %[0:1/fs:whole*(3/8)]; %dotted quarter note
ttrip = t1/12; %[0:1/fs:whole/12]; %triplet
t8 = t1/8; %[0:1/fs:whole/8]; %eighth note
t8dot = t1*3/16; %[0:1/fs:whole*(3/16)]; %dotted eighth note
t16 = t1/16; %[0:1/fs:whole/16]; %sixteenth note
t16dot = t1*5/16; %[0:1/fs:whole*(5/16)]; %dotted sixteenth note
rest8 = zeros(1,length(t8)); %eighth rest
rest16 = zeros(1,length(t16)); %sixteenth rest
bet = 0.001; %[0:1/fs:.001]; %length of pause b/w identical notes
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

b = zeros(1,bw*fs);

%%
m1 = [g(fs,Fsl,0.5) g(fs,A,0.5)];
m2 = [g(fs,B,t4) g(fs,A,t4) g(fs,Fsl,t4) zeros(1,t8*fs) g(fs,Fsl,t16) g(fs,A,t16)];
m3 = [g(fs,B,t4) g(fs,A,t4) g(fs,Fsl,t4) zeros(1,t8*fs) g(fs,Fsl,t16) g(fs,A,t16)];
m4 = [g(fs,B,t4) g(fs,Cs,t4) g(fs,D,t8) g(fs,A,t8dot) zeros(1,t16*fs) g(fs,Dl,t16) g(fs,Fsl,t16)];
m5 = [b g(fs,Fsl,t4) b g(fs,Fsl,ttrip) g(fs,Gl,ttrip) g(fs,Fsl,ttrip) g(fs,El,t4) zeros(1,t8*fs) g(fs,Fsl,t16) g(fs,A,t16)];
m6 = [g(fs,B,t4) g(fs,A,t16) b g(fs,A,t8) b g(fs,A,t16) g(fs,Fsl,t4dot) g(fs,D,t16) g(fs,Cs,t16)];
m7 = [g(fs,B,t4) g(fs,A,t8) b g(fs,A,t16) g(fs,Fsl,t16) g(fs,Dl,t4dot) g(fs,Fsl,t16) g(fs,A,t16)];
m8 = [g(fs,B,t4) g(fs,Cs,t4) g(fs,D,t16) g(fs,B,t8dot) g(fs,D,t16) g(fs,Fs,t8dot) b b];
m9 = [g(fs,Fs,t8) g(fs,G,t8) g(fs,Fs,t8) g(fs,G,t16) g(fs,E,t16) g(fs,E,t4) zeros(1,t8*fs) g(fs,Fs,t16) g(fs,G,t16)];
m10 = [0.8*g(fs,Ah,t4) g(fs,D,t8) g(fs,E,t2) zeros(1,t8*fs)];
m11 = [g(fs,Fs,t8) g(fs,E,t16) g(fs,Fs,t8dot) g(fs,G,t16) g(fs,Fs,t8), g(fs,E,t8dot) zeros(1,t8*fs) g(fs,Fs,t16) g(fs,G,t16)];
% %High A note is too high of a frequency, so reduced the amplitude
m12 = [0.8*g(fs,Ah,t4) g(fs,D,t8) g(fs,E,t2) zeros(1,t8*fs)];
m13 = [g(fs,Fs,t16) g(fs,E,t8) g(fs,Fs,t8) g(fs,G,t16) g(fs,Fs,t8) g(fs,E,t4) zeros(1,t8*fs) g(fs,Fs,t16) g(fs,G,t16)];
m14 = [0.8*g(fs,Ah,t8) g(fs,D,t8) g(fs,Cs,t8) g(fs,B,t2) g(fs,D,t16) g(fs,E,t16)];
m15 = [g(fs,Fs,t4) g(fs,A,t8) b g(fs,A,t2) zeros(1,t8*fs)];
m16 = [zeros(1,t8*fs) zeros(1,t8*fs) zeros(1,t8*fs) zeros(1,t8*fs) zeros(1,t8*fs) zeros(1,t8*fs) zeros(1,t8*fs) g(fs,Fsl,t16) g(fs,A,t16)];
m17 = [g(fs,B,t4) g(fs,A,t4) g(fs,Fsl,t8) g(fs,A,t8) g(fs,D,t8) g(fs,Fs,t8) b];
m18 = [g(fs,Fs,t4dot) g(fs,G,t8) g(fs,E,t8) g(fs,D,t16) g(fs,B,t8dot) zeros(1,t16*fs) g(fs,D,t16)];
m19 = [g(fs,D,t1)];

m = [m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13 m14 m15 m16 m17 m18 m19];
%soundsc(m,8000)
