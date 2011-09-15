% simulate periodic GP on a circle

close all
clear all

N = 100;
th = linspace(0,2*pi,N);

% th = th(1:end-1); % needed?

x = cos(th);
y = sin(th);
d = sqrt((x - x(1)).^2 + (y - y(1)).^2);

sigma = 1.0;
lambda = 0.25;
f = sigma^2*exp(-d.^2/(2*lambda^2));

% these are the eigenvalues
F = real(fft(f))/N;

% form the K-L expansion
Z = randn(N);
% N = 100: modes 1:51 are cosine, 2:50 are also sine
p = zeros(1,N);
for i = 1:N/2+1
    p = p + Z(i)*F(i)*cos(i*th);
end
for i = 2:N/2
    p = p + Z(i)*F(i)*sin(i*th);
end

plot(th,p)
hold on
plot(th+2*pi,p)
plot(th(end),p(end),'ro')

