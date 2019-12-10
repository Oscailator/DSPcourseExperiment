%%
clear
close all
clc
%%%%%%%%%Question 1, N=45���������δ����������������������Ĺ�һ��������%%%%%%%%%%%%%%%%%%
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
legend('���δ�', '������', '����������')
xlabel('t');ylabel('����');
title('���ִ��ں���');

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
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');
title('���ִ��ں���');
legend('���δ�', '������', '����������')

%%
clear
close all
clc
%%%%%%%%%Question 2, ��������ƴ�ͨ�˲���%%%%%%%%%%%%%%%%%%
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
legend('N=15,fir1���', 'N=45,fir1���')
xlim([-inf, pi])
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');

%%
clear
close all
clc
%%%%%%%%%Question 3, �����������δ���������������ƴ�ͨ�˲���%%%%%%%%%%%%%%%%%%
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
legend('������fir2' ,'���δ�fir2' ,'����������fir2', '������fir1' ,'���δ�fir1' ,'����������fir1')
title('N=15')
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');

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
legend('������fir1' ,'���δ�fir1' ,'����������fir1')
title('N=45')
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');

%%
clear
close all
clc
%%%%%%%%%Question 4, ���������������λ�˲���%%%%%%%%%%%%%%%%%%
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
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');
title('\beta=4�����������������λ�˲���');

phase=angle(h1);

subplot(3,2,2)
plot(w1/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('��һ��Ƶ��/\pi');
ylabel('��λ');

%%%%%%%%%%part 2:beta = 6%%%%%%%%
beta = 6;
h = fir2(N-1,f,a,kaiser(N,beta));
[h1,w1] = freqz(h,1);

subplot(3,2,3)
plot(w1/pi, 20*log10(abs(h1)), 'LineWidth' ,1.5);
%axis([0,1,-80,10]);
grid on
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');
title('\beta=6�����������������λ�˲���');

phase = angle(h1);

subplot(3,2,4)
plot(w1/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('��һ��Ƶ��/\pi');
ylabel('��λ');

%%%%%%%%%%part 3:beta = 10%%%%%%%%
beta = 10;
h = fir2(N-1,f,a,kaiser(N,beta));
[h1,w1] = freqz(h,1);
subplot(3,2,5)
plot(w1/pi, 20*log10(abs(h1)), 'LineWidth' ,1.5);
%axis([0,1,-120,10]);
grid on
xlabel('��һ��Ƶ��/\pi');
ylabel('����/dB');
title('\beta=10,���������������λ�˲���');

phase=angle(h1);

subplot(3,2,6)
plot(w1/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('��һ��Ƶ��/\pi');ylabel('��λ');

%%
clear
close all
clc
%%%%%%%%%Question 5, Ƶ�ʲ��������������λ�˲���%%%%%%%%%%%%%%%%%%
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
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');
title(' Ƶ�ʲ��������ר��������λ�˲���');

phase = angle(H);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('��һ��Ƶ��/\pi');ylabel('��λ');
title(' Ƶ�ʲ��������ר��������λ�˲���');

Hk = [zeros(1,3) 0.5 ones(1,5) 0.5 zeros(1,1) 0.5 ones(1,5) 0.5, zeros(1,5) 0.5 ones(1,5) 0.5 zeros(1,1) ones(1,5) 0.5 zeros(1,3)];
hn = imag(ifft(Hk.*exp(-1i*pi*(N-1)*k/N)));
[H, w] = freqz(hn,1);

figure()
subplot(2,1,1)
plot(w/pi, 20*log10(abs(H)), 'LineWidth' ,1.5);
axis([0 1 -80 10]);
grid on
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');
title(' Ƶ�ʲ��������ר��������λ�˲���');

phase = angle(H);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('��һ��Ƶ��/\pi');ylabel('��λ');
title(' Ƶ�ʲ��������ר��������λ�˲���');

%%
clear
close all
clc
%%%%%%%%%Question 6, �����ȷ����������λ�˲���%%%%%%%%%%%%%%%%%%
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
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');
title(' �����Ƚ����㷨ר��������λ�˲���');

phase = angle(h);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('��һ��Ƶ��/\pi');ylabel('��λ');
title(' �����Ƚ����㷨ר��������λ�˲���');

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
xlabel('��һ��Ƶ��/\pi');ylabel('����/dB');
title(' Ƶ�ʲ��������ר��������λ�˲���');

phase = angle(H);

subplot(2,1,2)
plot(w/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('��һ��Ƶ��/\pi');ylabel('��λ');

%%
clear
close all
clc
%%%%%%%%%Question 7, �����ȷ���Ƹ�ͨ������λ�˲���%%%%%%%%%%%%%%%%%%
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
xlabel('Ƶ��/\pi');
ylabel('����/dB');
title('�����Ƚ����㷨���������λ��ͨFIR�����˲���');

phase=angle(h);

subplot(2,1,2)
plot(w*2500/pi, phase, 'LineWidth' ,1.5);
grid on
xlabel('Ƶ��/\pi');ylabel('��λ');