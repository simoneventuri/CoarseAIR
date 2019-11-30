close all
clc
clear all

Tv = [2000, 3500, 5000, 7000];

C  = 2.5e14
n  = 0.0
Tr = 71000

K = C .* Tv.^n .* exp( - Tr ./ Tv) ./ 6.0221409e+23

figure
semilogy(10000./Tv,K)
hold on



C  = 2.5e14
n  = 0.0
Tr = 87740

K = C .* Tv.^n .* exp( - Tr ./ Tv) ./ 6.0221409e+23
semilogy(10000./Tv,K)




C  = 4.8e16
n  = -0.6
Tr = 40480

K = C .* Tv.^n .* exp( - Tr ./ Tv) ./ 6.0221409e+23
semilogy(10000./Tv,K)


