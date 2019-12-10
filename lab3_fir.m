%%
clear
close all
clc
%%%%%%%%%Question 1, N=45，画出矩形窗、汉明窗、布莱克曼窗的归一化幅度谱%%%%%%%%%%%%%%%%%%
N = 45;
n = 0:N-1;
W1 = boxcar(N);
W2 = hamming(N)';
W3 = blackman(N);

subplot(2,1,1)
plot(n, W1, 'LineWidth', 1.5)
hold on
plot(n, W2, '--', 'LineWidth',1.5)
hold on
plot(n, W3, '-.', 'LineWidth',1.5)
hold on
line([44, 44], [0, 1], 'LineWidth',1.5);
legend('矩形窗', '汉明窗', '布莱克曼窗')
xlabel('t');ylabel('幅度');
title('三种窗口函数');

[h1,w1]=freqz(W1,N);
[h2,w2]=freqz(W2,N);
[h3,w3]=freqz(W3,N);

subplot(2,1,2)
plot(w1/pi,20*log10(abs(h1)), 'LineWidth', 1.5)
hold on
plot(w2/pi,20*log10(abs(h2)), '-', 'LineWidth', 1.5)
hold on
plot(w3/pi,20*log10(abs(h3)), '-.', 'LineWidth', 1.5)
hold on
grid on
xlabel('归一化频率/\pi');ylabel('幅度/dB');
title('三种窗口函数');
legend('矩形窗', '汉明窗', '布莱克曼窗')

%%
clear
close all
clc
%%%%%%%%%Question 2, 汉宁窗设计带通滤波器%%%%%%%%%%%%%%%%%%
%%%%%%%%part 1: N = 15 %%%%%%%%%%%%%
N = 15;
ws = 2*pi;
wpl = 0.3*pi;
wph = 0.5*pi;

f = [wpl, wph];
f = f ./ ws * 2;
filterOut = fir1(N-1, f, 'bandpass', hanning(N));
[h15H, ~] = freqz(filterOut,1);
h15H = abs(h15H);
h15H = 20 * log10(h15H);

%%%%%%%%part 2: N = 45 %%%%%%%%%%%%%
N = 45;
f = [wpl, wph];
f = f ./ ws * 2;
filterOut = fir1(N-1, f, 'bandpass', hanning(N));
[h45, w] = freqz(filterOut,1);
h45 = abs(h45);
h45 = 20 * log10(h45);

%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%
plot(w, h15H, 'blue', 'LineWidth' ,1.5);
hold on
plot(w, h45, 'r', 'LineWidth' ,1.5);
grid on
legend('N=15,fir1设计', 'N=45,fir1设计')
xlim([-inf, pi])
xlabel('归一化频率/\pi');ylabel('幅度/dB');

%%
clear
close all
clc
%%%%%%%%%Question 3, 汉宁窗，矩形窗，布莱克曼窗设计带通滤波器%%%%%%%%%%%%%%%%%%
%%%%%%%%part 1: N = 15 %%%%%%%%%%%%%
N = 15;
ws = 2*pi;
wpl = 0.3*pi;
wph = 0.5*pi;

f = [wpl, wph];
f = f ./ ws * 2;
filterOut1 = fir1(N, f, 'bandpass', ones(1,N+1));
[h15H, ~] = freqz(filterOut1);
h15H = abs(h15H);
h15H = 20 * log10(h15H);

filterOut2 = fir1(N, f, 'bandpass', hanning(N+1));
[h15R, ~] = freqz(filterOut2);
h15R = abs(h15R);
h15R = 20 * log10(h15R);

filterOut3 = fir1(N, f, 'bandpass', blackman(N+1));
[h15B, w] = freqz(filterOut3);
h15B = abs(h15B);
h15B = 20 * log10(h15B);

subplot(2,1,1)
plot(w, h15H, 'blue', 'LineWidth' ,1.5);
hold on
plot(w, h15R, 'r', 'LineWidth' ,1.5);
hold on
plot(w, h15B, 'g', 'LineWidth' ,1.5);
grid on
legend('汉宁窗fir2' ,'矩形窗fir2' ,'布莱克曼窗fir2', '汉宁窗fir1' ,'矩形窗fir1' ,'布莱克曼窗fir1')
title('N=15')
xlabel('归一化频率/\pi');ylabel('幅度/dB');

%%%%%%%%part 2: N = 45 %%%%%%%%%%%%%
N = 45;

f = [wpl, wph];
f = f ./ ws * 2;
filterOut1 = fir1(N, f, 'bandpass', ones(1,N+1));
[h15H, ~] = freqz(filterOut1);
h15H = abs(h15H);
h15H = 20 * log10(h15H);

filterOut2 = fir1(N, f, 'bandpass', hanning(N+1));
[h15R, ~] = freqz(filterOut2);
h15R = abs(h15R);
h15R = 20 * log10(h15R);

filterOut3 = fir1(N, f, 'bandpass', blackman(N+1));
[h15B, w] = freqz(filterOut3);
h15B = abs(h15B);
h15B = 20 * log10(h15B);

subplot(2,1,2)
plot(w, h15H, 'blue', 'LineWidth' ,1.5);
hold on
plot(w, h15R, 'r', 'LineWidth' ,1.5);
hold on
plot(w, h15B, 'g', 'LineWidth' ,1.5);
grid on
legend('汉宁窗fir1' ,'矩形窗fir1' ,'布莱克曼窗fir1')
title('N=45')
xlabel('归一化频率/\pi');ylabel('幅度/dB');

%%
clear
close all
clc
%%%%%%%%%Question 4, 凯塞窗设计线性相位滤波器%%%%%%%%%%%%%%%%%%
N = 40;
f = [0 0.2 0.2 0.4 0.4 0.6 0.6 0.8 0.8 1];
a = [0 0 1 1 0 0 1 1 0 0];

%%%%%%%%%%part 1:beta = 4%%%%%%%%
beta = 4;
h = fir2(N-1,f,a,kaiser(N,beta));
[h1,w1] = freqz(h,1);

subplot(3,2,1)
plot(w1/pi, 20*log10(abs(h1)), 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');ylabel('幅度/dB');
title('\beta=4，凯塞窗设计线性相位滤波器');

phase=angle(h1);

subplot(3,2,2)
plot(w1/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');
ylabel('相位');

%%%%%%%%%%part 2:beta = 6%%%%%%%%
beta = 6;
h = fir2(N-1,f,a,kaiser(N,beta));
[h1,w1] = freqz(h,1);

subplot(3,2,3)
plot(w1/pi, 20*log10(abs(h1)), 'LineWidth' ,1.5);
%axis([0,1,-80,10]);
grid on
xlabel('归一化频率/\pi');ylabel('幅度/dB');
title('\beta=6，凯塞窗设计线性相位滤波器');

phase = angle(h1);

subplot(3,2,4)
plot(w1/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');
ylabel('相位');

%%%%%%%%%%part 3:beta = 10%%%%%%%%
beta = 10;
h = fir2(N-1,f,a,kaiser(N,beta));
[h1,w1] = freqz(h,1);
subplot(3,2,5)
plot(w1/pi, 20*log10(abs(h1)), 'LineWidth' ,1.5);
%axis([0,1,-120,10]);
grid on
xlabel('归一化频率/\pi');
ylabel('幅度/dB');
title('\beta=10,凯塞窗设计线性相位滤波器');

phase=angle(h1);

subplot(3,2,6)
plot(w1/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');ylabel('相位');

%%
clear
close all
clc
%%%%%%%%%Question 5, 频率采样法设计线性相位滤波器%%%%%%%%%%%%%%%%%%
N = 40;
k = 0:N-1;
Hk = [zeros(1,3) 0.5 ones(1,5) 0.5 zeros(1,1) 0.5 ones(1,5) 0.5, zeros(1,5) -0.5 -ones(1,5) -0.5 zeros(1,1) -ones(1,5) -0.5 zeros(1,3)];
hn = real(ifft(Hk.*exp(-1i*pi*(N-1)*k/N)));
[H, w] = freqz(hn,1);

figure(1)
subplot(2,1,1)
plot(w/pi, 20*log10(abs(H)), 'LineWidth' ,1.5);
axis([0 1 -80 10]);
grid on
xlabel('归一化频率/\pi');ylabel('幅度/dB');
title(' 频率采样法设计专用线性相位滤波器');

phase = angle(H);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');ylabel('相位');
title(' 频率采样法设计专用线性相位滤波器');

Hk = [zeros(1,3) 0.5 ones(1,5) 0.5 zeros(1,1) 0.5 ones(1,5) 0.5, zeros(1,5) 0.5 ones(1,5) 0.5 zeros(1,1) ones(1,5) 0.5 zeros(1,3)];
hn = imag(ifft(Hk.*exp(-1i*pi*(N-1)*k/N)));
[H, w] = freqz(hn,1);

figure()
subplot(2,1,1)
plot(w/pi, 20*log10(abs(H)), 'LineWidth' ,1.5);
axis([0 1 -80 10]);
grid on
xlabel('归一化频率/\pi');ylabel('幅度/dB');
title(' 频率采样法设计专用线性相位滤波器');

phase = angle(H);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');ylabel('相位');
title(' 频率采样法设计专用线性相位滤波器');

%%
clear
close all
clc
%%%%%%%%%Question 6, 雷米兹法设计线性相位滤波器%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%part 1%%%%%%%%%%%%%%
N = 40;
f = [0 0.15 0.2 0.4 0.45 0.55 0.6 0.8 0.85 1];
a = [0 0 1 1 0 0 1 1 0 0];
wt = [2 1 2 1 2];
b = remez(N-1,f,a,wt);
[h,w] = freqz(b,1);

figure(1)
subplot(2,1,1)
plot(w/pi, 20*log10(abs(h)), 'LineWidth' ,1.5);
% axis([0 1 -70 10])
grid on
xlabel('归一化频率/\pi');ylabel('幅度/dB');
title(' 雷米兹交替算法专用线性相位滤波器');

phase = angle(h);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');ylabel('相位');
title(' 雷米兹交替算法专用线性相位滤波器');

%%%%%%%%%%%%%%%%%%%part 2%%%%%%%%%%%%%%
k = 0:N-1;
Hk = [zeros(1,3) 0.5 ones(1,5) 0.5 zeros(1,1) 0.5 ones(1,5) 0.5, zeros(1,5) -0.5 -ones(1,5) -0.5 zeros(1,1) -ones(1,5) -0.5 zeros(1,3)];
hn = real(ifft(Hk.*exp(-1i*pi*(N-1)*k/N)));
[H, w] = freqz(hn,1);

figure(2)
subplot(2,1,1)
plot(w/pi, 20*log10(abs(H)), 'LineWidth' ,1.5);
axis([0 1 -80 10]);
grid on
xlabel('归一化频率/\pi');ylabel('幅度/dB');
title(' 频率采样法设计专用线性相位滤波器');

phase = angle(H);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('归一化频率/\pi');ylabel('相位');

%%
clear
close all
clc
%%%%%%%%%Question 7, 雷米兹法设计高通线性相位滤波器%%%%%%%%%%%%%%%%%%
fedge = [500 800];
mval = [0 1];
dev = [0.01 0.109];
fs = 5000;
[N, fpts, mag, wt] = remezord(fedge, mval, dev, fs);
b = remez(N, fpts, mag, wt);
[h, w] = freqz(b, 1);

subplot(2,1,1)
plot(w*2500/pi, 20*log10(abs(h)), 'LineWidth', 1.5);
axis([0 2500 -80 10])
grid on
xlabel('频率/\pi');
ylabel('幅度/dB');
title('雷米兹交替算法设计线性相位高通FIR数字滤波器');

phase=angle(h);

subplot(2,1,2)
plot(w*2500/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('频率/\pi');ylabel('相位');