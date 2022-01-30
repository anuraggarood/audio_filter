clc;
close all;
clear all;
%order=2;

%size=2;
%fs=8192;
%t=[0:1/fs:size];
%N=fs*size;
%% voice file
[y,f1] = audioread("Anurag G/audio/music.wav");
%TotalTime1 = length(y1)./fs1;
info = audioinfo("Anurag G/audio/music.wav");
voice = audioplayer(y,f1);
play(voice)
t = 0:seconds(1/f1):seconds(info.Duration);
t = t(1:end-1);
figure(1);
plot(t,y)
xlabel('Time')
ylabel('Audio Signal')
%% noise file
[x,f2] = audioread("Anurag G/audio/noise.wav");
info1 = audioinfo("Anurag G/audio/noise.wav");
noise = audioplayer(x,f2);
play(noise)
TotalTime2 = length(x)./f2;
t = (0:seconds(1/f2):seconds(TotalTime2));
t = t(1:end-1);
figure(2);
plot(t,x)
xlabel('Time')
ylabel('Noise Signal')
%% combine signal
%function [Y,f]=togau(music,noise)
[l1,c1]=size(y);
[l2,c2]=size(x);
if l1 > l2
    l=l2;
else
    l=l1;
end
if(f1>=f2)
    f=f1;
else
    f=f2;
end
G=zeros(1,1);
a=0;
token=0;
for i=1:l
    if a==f
        if token==0
            token=1;
            a=0;
        else
            token=0;
            a=0;
        end
    end
    if(token==0)
    G(i)=y(i);
    a=a+1;
    else
    G(i)=x(i);
    a=a+1;
    
    end
end
Y=G;
audiowrite('out.wav',Y,f);
figure(3);
subplot(2,1,1)
plot(t,Y)
xlabel('time')
ylabel('combine')
q = audioplayer(Y,f);
play(q)
%% combine echo and filtering it
num = [1,zeros(1,4800),0.9];
den = [1];
p = filter(num,den,y);
figure(4);
subplot(2,1,1)
plot(t,p)
xlabel('time')
ylabel('audio signal with echo')
echo = audioplayer(p,f)
play(echo)
%Removing echo
den=[1,zeros(1,4800),0.8];
num=[1];
r=filter(num,den,p);
z = audioplayer(r,f);
play(z);
TotalTime4 = length(r)./f;
t = (0:seconds(1/f):seconds(TotalTime4));
t = t(1:end-1);
subplot(2,1,2)
plot(t,r)
ylabel('audio(echo-filtered)')
xlabel('time')
%% combine siganl from audacity read and buffer filter
[y3,f3] = audioread('Anurag G/audio/combine.mp3');
combine = audioplayer(y3,f3);
play(combine)
TotalTime3 = length(y3)./f3;
t = (0:seconds(1/f3):seconds(TotalTime3));
t = t(1:end-1);
figure(5);
plot(t,y3)
xlabel('Time')
ylabel('combine Signal')
%% Plot both audio channels
N = size(y3,1); % Determine total number of samples in audio file
figure(6);
subplot(2,1,1);
stem(1:N, y3(:,1));
title('Left Channel');
subplot(2,1,2);
stem(1:N, y3(:,2));
title('Right Channel');
%% Plot the spectrum
df = f3 / N;
w = (-(N/2):(N/2)-1)*df;
y = fft(y3(:,1), N) / N; % For normalizing, but not needed for our analysis
y2 = fftshift(y);
figure(7);
plot(w,abs(y2));

% Design a bandpass filter that filters out between 700 to 12000 Hz
n = 7;
beginFreq = 700 / (f3/2);
endFreq = 12000 / (f3/2);
[b,a] = butter(n, [beginFreq, endFreq], 'bandpass');
%% Filter the signal
fOut = filter(b, a, y3);
TotalTime4 = length(fOut)./f3;
t = (0:seconds(1/f3):seconds(TotalTime4));
t = t(1:end-1);
figure(8);
plot(t,y3)
xlabel('Time')
ylabel('combine Signal (filtered)')
fil = audioplayer(fOut,f3);
play(fil)
%% filter design
audiowrite('echo-filtered.wav',r,f3)