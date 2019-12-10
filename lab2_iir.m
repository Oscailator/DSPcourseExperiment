%%
close all
clear
clc
%%%%%%%%%%%%%%�б�ѩ���ͨ�˲������%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 0.3;          % ͨ���߽�����0.3kHz
delta = 0.8;       % ͨ������0.8dB
fr = 0.2;          % ����߽�����0.2kHz
At = 20;           % ��С���˥��20dB
T = 1;             % ��������1ms
fs = 1/T;

% �б�ѩ��I��
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
%%%%%%%%%%%%%%������˹���ֵ�ͨ�˲������%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 0.2;          % ͨ���߽�����0.2kHz
delta = 1;         % ͨ������1dB
fr = 0.3;          % ����߽�����0.3kHz
At = 25;           % ��С���˥��25dB
T = 1;             % ��������1ms
fs = 1/T;

Wp = 2*pi*fc/fs;
Ws = 2*pi*fr/fs;

%%%%%%part 1: ������Ӧ���䷨%%%%%%%%%%
[N1, ~] = buttord(Wp, Ws, delta, At, 's');
[z1, p1 ,k1] = buttap(N1);
[B, A] = butter(N1, Wp, 's');
[num1, den1] = impinvar(B, A, 1);
[h1, w] = freqz(num1, den1);
% h1 = 20*log10(abs(h1));
f = w/(2*pi)*fs*1000;

%%%%%part 2: ˫���Ա任��%%%%%%%%%
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
legend(' ������Ӧ���䷨','˫���Ա任��')

figure(2)
zplane(p1,k1)
title('������Ӧ���䷨')

figure(3)
zplane(p2,k2)
title('˫���Ա任��')

%%
close all
clear
clc
%%%%%%%%%%%%%%˫���Ա仯����ư�����˹�ͣ��б�ѩ���ͺ���Բ�����ֵ�ͨ�˲������%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = 1.2;          % ͨ���߽�����1.2kHz
delta = 0.5;       % ͨ������0.5dB
fr = 2.0;          % ����߽�����2kHz
At = 40;           % ��С���˥��40dB
fs = 8.0;          % ����Ƶ��Ϊ8kHz
T = 1/fs;          % ��������(1/fs)ms

Wp = 2*pi*fc/fs;
Ws = 2*pi*fr/fs;
Wpb = 2/T*1000*tan(Wp/2);
Wsb = 2/T*1000*tan(Ws/2);

%%%%%%%%part 1:������˹��%%%%%%%%%%%%%%%%%
[N_butter, Wn_butter] = buttord(Wpb, Wsb, delta, At, 's');
[z_butter, p_butter ,k_butter] = buttap(N_butter);
[B, A] = butter(N_butter, Wpb, 's');
[num_butter, den_butter] = bilinear(B,A,1000*fs);
[h_butter, w_butter] = freqz(num_butter,den_butter);
f_butter = w_butter/(2*pi)*fs*1000;
h_butter = 20*log10(abs(h_butter));

%%%%%%%%part 2:�б�ѩ����%%%%%%%%%%%%%%%%%
[N_cheby1, Wn_cheby1] = cheb1ord(Wpb, Wsb, delta, At, 's');
[z_cheby1, p_cheby1, k_cheby1] = cheb1ap(N_cheby1, delta);
[B, A] = cheby1(N_cheby1, delta, Wn_cheby1, 'low', 's');
[num_cheby1, den_cheby1] = bilinear(B, A, 1000*fs);
[h_cheby1, w_cheby1] = freqz(num_cheby1,den_cheby1);
f_cheby1 = w_cheby1/(2*pi)*fs*1000;
h_cheby1 = 20*log10(abs(h_cheby1));

%%%%%%%%part 3:��Բ��%%%%%%%%%%%%%%%%%
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
legend(['������˹�ͣ�����Ϊ',num2str(N_butter)], ['�б�ѩ���ͣ�����Ϊ',num2str(N_cheby1)], ['��Բ�ͣ�����Ϊ',num2str(N_ellipap)])

figure(2)
zplane(p_butter,k_cheby1)
title('������˹��')

figure(3)
zplane(p_cheby1, k_ellip)
title('�б�ѩ����')

figure(4)
zplane(p_ellip, k_ellip)
title('��Բ��')
%%
close all
clear
clc
%%%%%%%%%%%%%%������˹���ִ�ͨ�˲������%%%%%%%%%%%%%%%%%%%%%%%%%%
fc_low = 2;          % ͨ���±߽�����2kHz
fc_high = 3;         % ͨ���ϱ߽�����3kHz
delta = 3;           % ͨ������3dB
fr_low = 1.5;        % ������߽�����1.5kHz
fr_high = 6;         % ������߽�����6kHz
At_low = 20;         % ��С�����˥��20dB
At_high = 5;         % ��С�����˥��20dB
fs = 30;             % ����Ƶ��Ϊ8kHz
T = 1/fs;            % ��������(1/fs)ms

Wp_low = 2*pi*fc_low/fs;
Wp_high = 2*pi*fc_high/fs;
Wp = [Wp_low,Wp_high];

Ws_low = 2*pi*fr_low/fs;
Ws_high = 2*pi*fr_high/fs;
Ws = [Ws_low,Ws_high];

At = max(At_low, At_high);


%%%%%%part 1: ������Ӧ���䷨%%%%%%%%%%
[N1, Wn1] = buttord(Wp, Ws, delta, At, 's');
[z1, p1 ,k1] = buttap(N1);
[B, A] = butter(N1, Wn1, 's');
[num1, den1] = impinvar(B, A, 1);
[h1, w] = freqz(num1, den1);
h1 = 20*log10(abs(h1));
f = w/(2*pi)*fs*1000;

%%%%%part 2: ˫���Ա任��%%%%%%%%%
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
legend(' ������Ӧ���䷨','˫���Ա任��')

figure(2)
zplane(p1, k1)
title('������Ӧ���䷨')

figure(3)
zplane(p2, k2)
title('˫���Ա任��')

%%
close all
clear
clc
%%%%%%%%%%%%%%˫���Ա仯������б�ѩ���ʹ����˲������%%%%%%%%%%%%%%%%%%%%%%%%%%
fc_low = 0.5;        % ͨ���±߽�����0.5kHz
fc_high = 3;         % ͨ���ϱ߽�����3kHz
delta = 3;           % ͨ������3dB
fr_low = 1;          % ������߽�����1kHz
fr_high = 2;         % ������߽�����2kHz
At = 18;             % ��С�����˥��20dB
fs = 10;             % ����Ƶ��Ϊ8kHz
T = 1/fs;            % ��������(1/fs)ms

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