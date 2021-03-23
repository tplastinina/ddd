function [retval] = main (z, r0, r1)
tic;
clc;
clear;
close all;
zm=350;
fm=350;
lambda=0.000532;
k=2*pi/lambda;
im=sqrt(-1); 
N=160;
r0=0.5;
r1=2.5;

x = linspace(-r1,r1,N);

y = x;

[x,y] = meshgrid(x,y);

r = sqrt(x.*x + y.*y);

fi = atan2(x,y);

fi_less_zero = find(fi<0);

fi(fi_less_zero) = fi(fi_less_zero) + 2*pi; 

Ak = zeros(size(x));

ringIdx = find((r>r0)&(r<r1));

Ak(ringIdx) = 1; 

Pg = zeros(size(x));

rnd = rand(1,length(ringIdx));

Pg(ringIdx) = rnd(1,:);  %cos(3.*fi)-2.*(sin(fi))^2;

g0 = Ak.*exp(1i*Pg);


%imshow(g0); %само ранд. кольцо
%figure;

G0=fft2(g0); %пф от рандомного кольца
aG0=abs(G0);
PG0=angle(G0); %фа3а ПФ от рандомного кольца

At1=double(imread ('bbb.png'));
At=At1(:,:,1);
size(At)
At=At/max(max(At)); %чтобы убрать большие числа при считывании
%At=At>=0.5; % к целому виду
%At=1-At; %белый на черном фоне

imshow(At); %треугольник
figure;

  
sG0=At.*exp(1i*PG0);  

sg0=ifft2(sG0); %обратное фурье
asg0=abs(sg0);
sPg0=angle(sg0);

g1=Ak.*exp(1i*sPg0);
  
G1=fft2(g1);
aG1=abs((G1));
PG1=angle((G1));

sG1 =At.*exp(1i*PG1);

sg1=ifft2(sG1);%обратное на 2 проходе
asg1=abs((sg1));
Psg1=angle((sg1));


g2=Ak .* exp(1i*Psg1);% на кольцо 


G2=fft2(g2);
aG2=abs(G2);
pG2=angle(G2);


sG2=At.*exp(1i*pG2);


sg2=ifft2(sG2); %обратное
asg2=abs(sg2);
psg2=angle(sg2);


g3 = Ak.*exp(1i*psg2);


G3=fft2(g3);
aG3=abs(G3);
pG3=angle(G3);

fri =exp(((1i*k)/2)*(x.*x + y.*y)*((1/zm)-(1/fm)));
%fri =exp(-(1i*k*(x.*x + y.*y))/(2*fm))*exp((1i*k*(x.*x + y.*y))/(2*zm)); %фринель множитель
frin = g1.*fri;
frinel=fft2(frin);
f=abs(frinel);
imagesc(f);%abs Фрeнель от g1
figure;

%imagesc(abs(G0)); %амплитуда пф от рандомного кольца
%figure;
%imagesc(PG0); %фаза пф от рандомного кольца
%figure;
%imagesc(asg0); %амплитуда обратного фурье от треугольника на е в ст.фа3ы рандомного кольца
%figure;
%imagesc(sPg0); %фаза обратного фурье от тр. на е в ст. фа3ы рандомного кольца
%figure;
imagesc(aG1);
%figure;
%imagesc(PG1);
%figure();
f_amp = max(f(:));
f(f <= (f_amp*0.5)) = 0; 

[imgDiff,imgDiffAvg,imgDiffRms] = ImgSimilarity(f,At);
imagesc(imgDiff);
imgDiffRms

%imagesc(asg1);
%figure;
%imagesc(Psg1);
%figure;
%imagesc(aG2);
%figure;
%imagesc(pG2);
%figure;
%imagesc(asg2);
%figure;
%imagesc(psg2);
%figure;
%imagesc(aG3);
%figure;
%imagesc(pG3);

%s=aG3-aG0;
%imagesc(s);
toc;
  end
