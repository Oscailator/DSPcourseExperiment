%%
close all
clear
clc
%%%%%%%%%%%%%%切比雪夫高通滤波器设计%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 0.3;          % 通带边界条件0.3kHz
delta = 0.8;       % 通带波动0.8dB
fr = 0.2;          % 阻带边界条件0.2kHz
At = 20;           % 最小阻带衰减20dB
T = 1;             % 采样周期1ms
fs = 1/T;

% 切比雪夫I型
Wp = 2*pi*fc/fs;
Ws = 2*pi*fr/fs;
[N, Wn] = cheb1ord(Wp,Ws, delta, At, 's');
[z, p, k] = cheb1ap(N, delta);
[B, A] = cheby1(N, delta, Wn, 'high', 's');
omega = 0:0.01:2*pi;
h = freqz(B, A, omega);
gain = 20*log10(abs(h));

figure(1)
plot(omega/(2*pi),gain,'LineWidth',1.5)
ylim([-20,5])
grid on

figure(2)
zplane(p,k)

%%
close all
clear
clc
%%%%%%%%%%%%%%巴特沃斯数字低通滤波器设计%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 0.2;          % 通带边界条件0.2kHz
delta = 1;         % 通带波动1dB
fr = 0.3;          % 阻带边界条件0.3kHz
At = 25;           % 最小阻带衰减25dB
T = 1;             % 采样周期1ms
fs = 1/T;

Wp = 2*pi*fc/fs;
Ws = 2*pi*fr/fs;

%%%%%%part 1: 脉冲响应不变法%%%%%%%%%%
[N1, ~] = buttord(Wp, Ws, delta, At, 's');
[z1, p1 ,k1] = buttap(N1);
[B, A] = butter(N1, Wp, 's');
[num1, den1] = impinvar(B, A, 1);
[h1, w] = freqz(num1, den1);
% h1 = 20*log10(abs(h1));
f = w/(2*pi)*fs*1000;

%%%%%part 2: 双线性变换法%%%%%%%%%
Wpb = 2/T*1000*tan(Wp/2);
Wsb = 2/T*1000*tan(Ws/2);
[N2, ~] = buttord(Wpb, Wsb, delta, At, 's');
[z2, p2 ,k2] = buttap(N2);
[B, A] = butter(N2, Wpb, 's');
[num2, den2] = bilinear(B,A,1000*fs);
[h2, ~] = freqz(num2,den2);

%%%%plot%%%%%
figure(1)
plot(f, abs(h1), '-.', 'LineWidth', 1.5)
hold on
plot(f, abs(h2), 'LineWidth', 1.5)
grid on
legend(' 脉冲响应不变法','双线性变换法')

figure(2)
zplane(p1,k1)
title('脉冲响应不变法')

figure(3)
zplane(p2,k2)
title('双线性变换法')

%%
close all
clear
clc
%%%%%%%%%%%%%%双线性变化法设计巴特沃斯型，切比雪夫型和椭圆型数字低通滤波器设计%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 1.2;          % 通带边界条件1.2kHz
delta = 0.5;       % 通带波动0.5dB
fr = 2.0;          % 阻带边界条件2kHz
At = 40;           % 最小阻带衰减40dB
fs = 8.0;          % 采样频率为8kHz
T = 1/fs;          % 采样周期(1/fs)ms

Wp = 2*pi*fc/fs;
Ws = 2*pi*fr/fs;
Wpb = 2/T*1000*tan(Wp/2);
Wsb = 2/T*1000*tan(Ws/2);

%%%%%%%%part 1:巴特沃斯型%%%%%%%%%%%%%%%%%
[N_butter, Wn_butter] = buttord(Wpb, Wsb, delta, At, 's');
[z_butter, p_butter ,k_butter] = buttap(N_butter);
[B, A] = butter(N_butter, Wpb, 's');
[num_butter, den_butter] = bilinear(B,A,1000*fs);
[h_butter, w_butter] = freqz(num_butter,den_butter);
f_butter = w_butter/(2*pi)*fs*1000;
h_butter = 20*log10(abs(h_butter));

%%%%%%%%part 2:切比雪夫型%%%%%%%%%%%%%%%%%
[N_cheby1, Wn_cheby1] = cheb1ord(Wpb, Wsb, delta, At, 's');
[z_cheby1, p_cheby1, k_cheby1] = cheb1ap(N_cheby1, delta);
[B, A] = cheby1(N_cheby1, delta, Wn_cheby1, 'low', 's');
[num_cheby1, den_cheby1] = bilinear(B, A, 1000*fs);
[h_cheby1, w_cheby1] = freqz(num_cheby1,den_cheby1);
f_cheby1 = w_cheby1/(2*pi)*fs*1000;
h_cheby1 = 20*log10(abs(h_cheby1));

%%%%%%%%part 3:椭圆型%%%%%%%%%%%%%%%%%
[N_ellipap, Wn_ellip] = ellipord(Wpb, Wsb, delta, At, 's');
[z_ellip, p_ellip, k_ellip] = ellipap(N_ellipap, delta, At);
[B, A] = ellip(N_ellipap, delta, At, Wn_ellip, 's');
[num_ellipap, den_ellipap] = bilinear(B, A, 1000*fs);
[h_ellipap, w_ellipap] = freqz(num_ellipap,den_ellipap);
f_ellipap = w_ellipap/(2*pi)*fs*1000;
h_ellipap = 20*log10(abs(h_ellipap));

%%%%%%plot%%%%%%%%
figure(1)
plot(f_butter, h_butter, '-.', 'LineWidth', 1.5)
hold on
plot(f_cheby1, h_cheby1, 'LineWidth', 1.5)
hold on
plot(f_ellipap, h_ellipap, '--', 'LineWidth', 1.5)
grid on
ylim([-50,1])
legend(['巴特沃斯型，阶数为',num2str(N_butter)], ['切比雪夫型，阶数为',num2str(N_cheby1)], ['椭圆型，阶数为',num2str(N_ellipap)])

figure(2)
zplane(p_butter,k_cheby1)
title('巴特沃斯型')

figure(3)
zplane(p_cheby1, k_ellip)
title('切比雪夫型')

figure(4)
zplane(p_ellip, k_ellip)
title('椭圆型')
%%
close all
clear
clc
%%%%%%%%%%%%%%巴特沃斯数字带通滤波器设计%%%%%%%%%%%%%%%%%%%%%%%%%%
fc_low = 2;          % 通带下边界条件2kHz
fc_high = 3;         % 通带上边界条件3kHz
delta = 3;           % 通带波动3dB
fr_low = 1.5;        % 下阻带边界条件1.5kHz
fr_high = 6;         % 上阻带边界条件6kHz
At_low = 20;         % 最小下阻带衰减20dB
At_high = 5;         % 最小上阻带衰减20dB
fs = 30;             % 采样频率为8kHz
T = 1/fs;            % 采样周期(1/fs)ms

Wp_low = 2*pi*fc_low/fs;
Wp_high = 2*pi*fc_high/fs;
Wp = [Wp_low,Wp_high];

Ws_low = 2*pi*fr_low/fs;
Ws_high = 2*pi*fr_high/fs;
Ws = [Ws_low,Ws_high];

At = max(At_low, At_high);


%%%%%%part 1: 脉冲响应不变法%%%%%%%%%%
[N1, Wn1] = buttord(Wp, Ws, delta, At, 's');
[z1, p1 ,k1] = buttap(N1);
[B, A] = butter(N1, Wn1, 's');
[num1, den1] = impinvar(B, A, 1);
[h1, w] = freqz(num1, den1);
h1 = 20*log10(abs(h1));
f = w/(2*pi)*fs*1000;

%%%%%part 2: 双线性变换法%%%%%%%%%
Wpb = 2/T*1000*tan(Wp/2);
Wsb = 2/T*1000*tan(Ws/2);

[N2, Wn2] = buttord(Wpb, Wsb, delta, At, 's');
[z2, p2 ,k2] = buttap(N2);
[B, A] = butter(N2, Wn2, 's');
[num2, den2] = bilinear(B,A,1000*fs);
[h2, ~] = freqz(num2,den2);
h2 = 20*log10(abs(h2));

%%%%plot%%%%%
figure(1)
plot(f, h1, '-.', 'LineWidth', 1.5)
hold on
plot(f, h2, 'LineWidth', 1.5)
grid on
xlim([0,15000])
ylim([-100, 10])
legend(' 脉冲响应不变法','双线性变换法')

figure(2)
zplane(p1, k1)
title('脉冲响应不变法')

figure(3)
zplane(p2, k2)
title('双线性变换法')

%%
close all
clear
clc
%%%%%%%%%%%%%%双线性变化法设计切比雪夫型带阻滤波器设计%%%%%%%%%%%%%%%%%%%%%%%%%%
fc_low = 0.5;        % 通带下边界条件0.5kHz
fc_high = 3;         % 通带上边界条件3kHz
delta = 3;           % 通带波动3dB
fr_low = 1;          % 下阻带边界条件1kHz
fr_high = 2;         % 上阻带边界条件2kHz
At = 18;             % 最小下阻带衰减20dB
fs = 10;             % 采样频率为8kHz
T = 1/fs;            % 采样周期(1/fs)ms

Wp_low = 2*pi*fc_low/fs;
Wp_high = 2*pi*fc_high/fs;
Wp = [Wp_low,Wp_high];

Ws_low = 2*pi*fr_low/fs;
Ws_high = 2*pi*fr_high/fs;
Ws = [Ws_low,Ws_high];

Wpb = 2/T*1000*tan(Wp/2);
Wsb = 2/T*1000*tan(Ws/2);

[N, Wn] = cheb1ord(Wpb, Wsb, delta, At, 's');
[z, p, k] = cheb1ap(N, delta);
[B, A] = cheby1(N, delta, Wn, 'stop', 's');
[num, den] = bilinear(B, A, 1000*fs);
[h, w] = freqz(num,den);
f = w/(2*pi)*fs*1000;
h = 20*log10(abs(h));

figure(1)
plot(f, h, 'LineWidth', 1.5)
grid on

figure(2)
zplane(p, k)