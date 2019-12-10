%%
close all
clear
clc
%%%%%%%%%%%%%%%question 1:观察高斯序列的时域和幅频特性%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Part 1 fixed p%%%%%%%%%%%%%
n = 0:15;
p = 8;
fig_num = 1;
for q = [2,4,8]
    x = GuassionSeq(n, p, q);
    y = abs(fft(x));
    
    figure(fig_num)
    subplot(1,2,1)
    stem(n, x, 'LineWidth', 1.5)
    title(['p=',num2str(p),' q=',num2str(q),' 时域图'])
    grid on
    subplot(1,2,2)
    stem(0:(length(y)-1), y, 'LineWidth', 1.5)
    title(['p=',num2str(p),' q=',num2str(q),' 幅频图'])
    grid on
    fig_num = fig_num + 1;
end

%%%%%%Part 2 fixed q%%%%%%%%%%%%%
q = 8;
for p = [8,13,14]
    x = GuassionSeq(n, p, q);
    y = abs(fft(x));   
    figure(fig_num)
    subplot(1,2,1)
    stem(n, x, 'LineWidth', 1.5)
    title(['p=',num2str(p),'q=',num2str(q),'时域图'])
    grid on
    subplot(1,2,2)
    stem(0:(length(y)-1), y, 'LineWidth', 1.5)
    title(['p=',num2str(p),'q=',num2str(q),'幅频图'])
    grid on
    fig_num = fig_num + 1;
end

%%
close all
clear
clc
%%%%%%%%%%%%%%%question 2:观察衰减正弦序列的时域和幅频特性%%%%%%%%%%%%%%%%%%%%%%
n = 0:15;
a = 0.1;
fig_num = 1;
for f = [0.0625,0.4375,0.5625]
    x = attenuationSin(n, a, f);
    y = abs(fft(x));
    
    figure(fig_num)
    subplot(1,2,1)
    stem(n, x, 'LineWidth', 1.5)
    title(['a=',num2str(a),' f=',num2str(f),' 时域图'])
    grid on
    subplot(1,2,2)
    stem(0:(length(y)-1), y, 'LineWidth', 1.5)
    title(['a=',num2str(a),' f=',num2str(f),' 幅频图'])
    grid on
    fig_num = fig_num + 1;
end

%%
close all
clear
clc
%%%%%%%%%%%%%%%question 3:观察三角序列和反三角序列的时域和幅频特性%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Part 1 三角序列%%%%%%%%%%%%%
n = 0:7;
x1 = Triangle(n);
y1 = abs(fft(x1));

figure(1)
subplot(2,2,1)
stem(n, x1, 'LineWidth', 1.5)
grid on
title('三角时域图')
subplot(2,2,2)
stem(0:(length(y1)-1), y1, 'LineWidth', 1.5)
grid on
title('三角幅频图')

%%%%%%Part 2 反三角序列%%%%%%%%%%%%%
x2 = reverseTriangle(n);
y2 = abs(fft(x2));

subplot(2,2,3)
stem(n, x2, 'LineWidth', 1.5)
grid on
title('反三角时域图')
subplot(2,2,4)
stem(0:(length(y2)-1), y2, 'LineWidth', 1.5)
grid on
title('反三角幅频图')

%%%%%%Part 3 补零三角序列%%%%%%%%%%%%%
n = 0:31;
x1 = Triangle(n);
y1 = abs(fft(x1));

figure(2)
subplot(2,2,1)
stem(n, x1, 'LineWidth', 1.5)
grid on
title('三角时域图')
subplot(2,2,2)
stem(0:(length(y1)-1), y1, 'LineWidth', 1.5)
grid on
title('三角幅频图')

%%%%%%Part 4 补零反三角序列%%%%%%%%%%%%%
x2 = reverseTriangle(n);
y2 = abs(fft(x2));

subplot(2,2,3)
stem(n, x2, 'LineWidth', 1.5)
grid on
title('反三角时域图')
subplot(2,2,4)
stem(0:(length(y2)-1), y2, 'LineWidth', 1.5)
grid on
title('反三角幅频图')

%%
close all
clear
clc
%%%%%%%%%%%%%%%question 4:观察信号长度对fft的影响%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%part 1: N = 16%%%%%%%
n = 0:15;
fignum = 1;
for deltaf = [1/16,1/64]
    x = sin(2*pi*0.125*n) + sin(2*pi*(0.125+deltaf)*n);
    y = abs(fft(x));
    
    figure(fignum)
    subplot(2,1,1)
    stem(n, x, 'LineWidth', 1.5)
    title(['16点 频率相差为',num2str(deltaf),'时域图'])
    grid on
    
    subplot(2,1,2)
    stem(0:(length(y)-1), y, 'LineWidth', 1.5)
    title(['16点 频率相差为',num2str(deltaf),'幅频图'])
    grid on
    fignum = fignum + 1;
end

%%%%%%%part 2: N = 128%%%%%%%
n = 0:127;
for deltaf = [1/16,1/64]
    x = sin(2*pi*0.125*n) + sin(2*pi*(0.125+deltaf)*n);
    y = abs(fft(x));
    
    figure(fignum)
    subplot(2,1,1)
    stem(n, x, 'LineWidth', 1.5)
    title(['128点 频率相差为',num2str(deltaf),'时域图'])
    grid on
    
    subplot(2,1,2)
    stem(0:(length(y)-1), y, 'LineWidth', 1.5)
    title(['128点 频率相差为',num2str(deltaf),'幅频图'])
    grid on
    fignum = fignum + 1;
end

%%
close all
clear
clc
%%%%%%%%%%%%%%%question 5:FFT计算卷积%%%%%%%%%%%%%%%%%%%%%%
seqLen = 16;
n = 0:seqLen-1;
xa = GuassionSeq(n, 8, 2);
xb = attenuationSin(n, 0.1, 0.0625);

figure(1)
stem(n, xa, 'LineWidth', 1.5)
grid on
title('高斯序列时域图')

figure(2)
stem(n, xb, 'LineWidth', 1.5)
grid on
title('衰减三角序列时域图')
%%%%%%%part 1: 循环卷积%%%%%%%
y1 = ifft(fft(xa).*fft(xb));
figure(3)
subplot(3,1,1)
stem(0:(length(y1)-1), y1, 'LineWidth', 1.5)
grid on
title('循环卷积时域图')
%%%%%%%part 2: 线性卷积%%%%%%%
y2 = ifft(fft(xa, 2*seqLen-1).*fft(xb, 2*seqLen-1));
subplot(3,1,2)
stem(0:(length(y2)-1), y2, 'LineWidth', 1.5)
grid on
title('线性卷积时域图')
%%%%%%part 3:fft自带卷积%%%%%%
y3 = conv(xa, xb);
subplot(3,1,3)
stem(0:(length(y3)-1), y3, 'LineWidth', 1.5)
grid on
title('matlab内置线性卷积时域图')

%%
close all
clear
clc
%%%%%%%%%%%%%%%question 6:随机序列卷积三角序列%%%%%%%%%%%%%%%%%%%%%%
N = 8;
xc = Triangle(0:7);
M = 512;
xe = randn(1,512);
y_conv = conv(xc,xe);

figure(1)
subplot(2,1,1)
stem(0:511, xe, 'LineWidth', 1.5)
grid on
title('随机序列时域图')

fbefore = abs(fft(xe));
subplot(2,1,2)
stem(0:length(fbefore)-1, fbefore, 'LineWidth', 1.5)
grid on
title('卷积前频谱')
%%%%%%%重叠相加法%%%%%%
numP = 8;
NP = M/numP;
L = N + NP - 1;
H = fft(xc,L);
resCutted1 = zeros(1,N+M-1);
for ii = 1:NP:M
    temp = ifft(fft(xe(1,ii:ii+NP-1), L) .* H);
    resCutted1(1,ii:ii+L-1) = resCutted1(1,ii:ii+L-1) + temp;
end

figure(2)
subplot(2,2,1)
stem(0:length(resCutted1)-1, resCutted1, 'LineWidth', 1.5)
grid on
title('重叠相加法卷积后时域图')

fatfer1 = abs(fft(resCutted1));
subplot(2,2,2)
stem(0:length(fatfer1)-1, fatfer1, 'LineWidth', 1.5)
grid on
title('重叠相加法卷积后频谱')

if(max(y_conv-resCutted1) > 1e-10)
    disp('wrong');
else
    disp('right');
end

%%%%%%%重叠保留法%%%%%%
resCutted2 = zeros(1,N+M-1);
H = fft(xc,L);
for ii = 1:NP:M
    if ii < N - 1
        temp = [zeros(1,N-1), xe(1:NP)];
    else
        temp = xe((ii-N+1):(ii+NP-1));
    end
    temp = ifft(fft(temp,L) .* H);
    resCutted2(1,ii:ii+NP-1) = temp(N:L);
end
temp = [xe(M-N+1:M),zeros(1,L-N+1)];
temp = ifft(fft(temp,L) .* H);
resCutted2(1,M:N+M-1) = temp(N:N+N-1);

subplot(2,2,3)
stem(0:length(resCutted2)-1, resCutted2, 'LineWidth', 1.5)
grid on
title('重叠保留法卷积后时域图')

fatfer2 = abs(fft(resCutted2));
subplot(2,2,4)
stem(0:length(fatfer2)-1, fatfer2, 'LineWidth', 1.5)
grid on
title('重叠保留法卷积后频谱')

if(max(y_conv-resCutted2) > 1e-10)
    disp('wrong');
else
    disp('right');
end

%%
close all
clear
clc
%%%%%%%%%%%%%%%question 7:FFT计算互相关%%%%%%%%%%%%%%%%%%%%%%
seqLen = 16;
n = 0:seqLen-1;
xa = GuassionSeq(n, 8, 2);
xb = attenuationSin(n, 0.1, 0.0625);

figure(1)
stem(n, xa, 'LineWidth', 1.5)
grid on
title('高斯序列时域图')

figure(2)
stem(n, xb, 'LineWidth', 1.5)
grid on
title('衰减三角序列时域图')
%%%%%%%part 1: 循环相关%%%%%%%
y1 = ifft(conj(fft(xa)).*fft(xb));
y1 = [y1(seqLen/2+2:seqLen), y1(1:seqLen/2)];
figure(3)
subplot(3,1,1)
stem(0:(length(y1)-1), y1(1:length(y1)), 'LineWidth', 1.5)
grid on
title('Rba循环相关图')
%%%%%%%part 2: 线性相关%%%%%%%
y2 = ifft(conj(fft(xa, 2*seqLen)).*fft(xb, 2*seqLen));
y2 = [y2(seqLen+2:2*seqLen), y2(1:seqLen)];
subplot(3,1,2)
stem(0:(length(y2)-1), y2, 'LineWidth', 1.5)
grid on
title('Rba线性相关图')
%%%%%%%part 2: 线性相关%%%%%%%
y4 = xcorr(xb, xa);
subplot(3,1,3)
stem(0:(length(y4)-1), y4, 'LineWidth', 1.5)
grid on
title('matlab内置Rba线性相关图')

%%%%%交换ab顺序%%%%%%%
%%%%%%%part 1: 循环相关%%%%%%%
y1 = ifft(conj(fft(xb)).*fft(xa));
% y1 = fftshift(y1);
y1 = [y1(seqLen/2+2:seqLen), y1(1:seqLen/2)];
figure(4)
subplot(3,1,1)
stem(0:(length(y1)-1), y1(1:(length(y1))), 'LineWidth', 1.5)
grid on
title('Rab循环相关图')
%%%%%%%part 2: 线性相关%%%%%%%
y2 = ifft(conj(fft(xb, 2*seqLen)).*fft(xa, 2*seqLen));
y2 = [y2(seqLen+2:2*seqLen), y2(1:seqLen)];
subplot(3,1,2)
stem(0:(length(y2)-1), y2, 'LineWidth', 1.5)
grid on
title('Rab线性相关图')
%%%%%%%part 2: 线性相关%%%%%%%
y4 = xcorr(xa, xb);
subplot(3,1,3)
stem(0:(length(y4)-1), y4, 'LineWidth', 1.5)
grid on
title('matlab内置Rab线性相关图')

%%
close all
clear
clc
%%%%%%%%%%%%%%%question 8:FFT计算自相关%%%%%%%%%%%%%%%%%%%%%%
seqLen = 16;
n = 0:seqLen-1;
xa = GuassionSeq(n, 8, 2);
xb = attenuationSin(n, 0.1, 0.0625);

figure(1)
subplot(2,2,1)
stem(n, xa, 'LineWidth', 1.5)
grid on
title('高斯序列时域图')

subplot(2,2,2)
stem(n, xb, 'LineWidth', 1.5)
grid on
title('衰减三角序列时域图')
%%%%%%%part 1: 高斯序列%%%%%%%%%%%
temp = fft(xa, 2*seqLen);
y1 = ifft(conj(temp).*temp);
y1 = fftshift(y1);

subplot(2,2,3)
stem(0:(length(y1)-1), y1, 'LineWidth', 1.5)
grid on
title('高斯序列自相关图')
%%%%%%%part 2: 衰减三角序列%%%%%%%
temp = fft(xb, 2*seqLen);
y2 = ifft(conj(temp) .* temp);
y3 = fftshift(y2);
y2 = [y2(seqLen+2:2*seqLen), y2(1:seqLen)];

subplot(2,2,4)
stem(0:(length(y2)-1), y2, 'LineWidth', 1.5)
grid on
title('衰减三角序列自相关图')