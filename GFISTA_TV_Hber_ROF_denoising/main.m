clear all
close all
clc
format longe
%%
% iter1=input('iteration1 = ');
% iter2=input('iteration2 = ');
iter1=5000;
iter2=100;
L_f=8;
mu=0.1;
mu_g=0.1;
mu_f=0;
lambda=0.1;
ep=0.01;
rho=0.9;
L_0_5=5;
L_0_20=20;
%%
u = im2double(imread('peppers256_est_s5.png'));
[m, n] = size(u);
% Add Gaussian noise
sigma = sqrt(0.005);
noise = sigma * randn(m, n);
u0 = u + noise;
% u0=imnoise(img,'gaussian',0,0.005);  %or
%%
[ u_star, F_u, improvement ]=GFISTA( u, u0, iter1, L_f, mu, mu_f, mu_g, lambda, ep );
Fu_star=F_u(iter1+1);
[ u1_1, L1_1, e1_1 ]=GFISTA_modified( u0, Fu_star, iter2, L_f, mu, mu_f, mu_g, lambda, ep );
[ u1_2_5, L1_2_5, e1_2_5 ]=GFISTA_cl_backtracking( u0, Fu_star, iter2, mu, mu_f, mu_g, lambda, ep, L_0_5, rho );
[ u1_3_5, L1_3_5, e1_3_5 ]=GFISTA_full_backtracking( u0, Fu_star, iter2, mu, mu_f, mu_g, lambda, ep, L_0_5, rho );
[ u1_4, e1_4 ] = FISTA ( u0, Fu_star, iter2, L_f, lambda, ep );
[ u1_2_20, L1_2_20, e1_2_20 ]=GFISTA_cl_backtracking( u0, Fu_star, iter2, mu, mu_f, mu_g, lambda, ep, L_0_20, rho );
[ u1_3_20, L1_3_20, e1_3_20 ]=GFISTA_full_backtracking( u0, Fu_star, iter2, mu, mu_f, mu_g, lambda, ep, L_0_20, rho );
%%
disp('------------------------------')
disp(['iteration = ', num2str(iter2)])
disp(['L_f = 8  >> Error_GFISTA = ',num2str(e1_1(iter2+1))])
disp(['L_0 = 5  >> Error_GFISTA_cl_bactracking = ',num2str(e1_2_5(iter2+1))])
disp(['L_0 = 5  >> Error_GFISTA_full_backtracking = ',num2str(e1_3_5(iter2+1))])
disp(['L_f = 8  >> Error_FISTA = ',num2str(e1_4(iter2+1))])
disp(['L_0 = 20 >> Error_GFISTA_cl_backtracking = ',num2str(e1_2_20(iter2+1))])
disp(['L_0 = 20 >> Error_GFISTA_full_backtracking = ',num2str(e1_3_20(iter2+1))])
%%
figure(1);
%-----
subplot(3,3,1);
imshow(u);
title('Original image.');
%-----
subplot(3,3,2);
imshow(u0);
title('Noisy version.');
%-----
subplot(3,3,3);
imshow(u_star);
title(['Denoised, GFISTA, iter=', num2str(iter1)]);
%-----
subplot(3,3,4);
imshow(u1_1);
title('GFISTA modified');
%-----
subplot(3,3,5);
imshow(u1_2_5);
title('GFISTA cl backtracking L_0=5');
%-----
subplot(3,3,6);
imshow(u1_3_5);
title('GFISTA full backtracking L_0=5');
%-----
subplot(3,3,7);
imshow(u1_4);
title('FISTA ');
%-----
subplot(3,3,8);
imshow(u1_2_20);
title('GFISTA cl backtracking L_0=20');
%-----
subplot(3,3,9);
imshow(u1_3_20);
title('GFISTA full backtracking L_0=20');
%%
figure(2);
subplot(2,3,1);
k=0:(round(iter1/5));
plot(k,F_u(k+1),'-b');
xlabel('k')
ylabel('F(u)')
%-------------------
subplot(2,3,2);
k=(round(iter1/5)):((2*round(iter1/5)));
plot(k,F_u(k+1),'-b');
xlabel('k')
ylabel('F(u)')
%-------------------
subplot(2,3,3);
k=((2*round(iter1/5))):((3*round(iter1/5)));
plot(k,F_u(k+1),'-b');
xlabel('k')
ylabel('F(u)')
%-------------------
subplot(2,3,4);
k=((3*round(iter1/5))):((4*round(iter1/5)));
plot(k,F_u(k+1),'-b');
xlabel('k')
ylabel('F(u)')
%-------------------
subplot(2,3,5);
k=((4*round(iter1/5))):((5*round(iter1/5)));
plot(k,F_u(k+1),'-b');
xlabel('k')
ylabel('F(u)')
%-------------------
subplot(2,3,6);
k=0:iter1;
plot(k,F_u(k+1),'-b');
xlabel('k')
ylabel('F(u)')
%%
figure(3)
k=1:(iter1+1);
plot(k,improvement(k),'-b');
title('improvement');
xlabel('k')
ylabel('improvement(k)')
%%
subplot(2,3,1);
k=1:(round(iter1/5));
plot(k,improvement(k),'-b');
xlabel('k')
ylabel('improvement(k)')
%-------------------
subplot(2,3,2);
k=(round(iter1/5)):((2*round(iter1/5)));
plot(k,improvement(k),'-b');
xlabel('k')
ylabel('improvement(k)')
%-------------------
subplot(2,3,3);
k=((2*round(iter1/5))):((3*round(iter1/5)));
plot(k,improvement(k),'-b');
xlabel('k')
ylabel('improvement(k)')
%-------------------
subplot(2,3,4);
k=((3*round(iter1/5))):((4*round(iter1/5)));
plot(k,improvement(k),'-b');
xlabel('k')
ylabel('improvement(k)')
%-------------------
subplot(2,3,5);
k=((4*round(iter1/5))):((5*round(iter1/5)));
plot(k,improvement(k),'-b');
xlabel('k')
ylabel('improvement(k)')
%-------------------
subplot(2,3,6);
k=1:iter1;
plot(k,improvement(k),'-b');
xlabel('k')
ylabel('improvement(k)')
%%
figure(4)
subplot(2,3,1);
k=0:20;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_5(k+1),'-r');
hold on
plot(k,e1_3_5(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=5');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,2);
k=20:40;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_5(k+1),'-r');
hold on
plot(k,e1_3_5(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=5');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,3);
k=40:60;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_5(k+1),'-r');
hold on
plot(k,e1_3_5(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=5');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,4);
k=60:80;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_5(k+1),'-r');
hold on
plot(k,e1_3_5(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=5');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,5);
k=80:iter2;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_5(k+1),'-r');
hold on
plot(k,e1_3_5(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=5');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%------------------------------------------------
subplot(2,3,6);
k=0:iter2;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_5(k+1),'-r');
hold on
plot(k,e1_3_5(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=5');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%%
figure(5)
k=0:iter2;
plot(k,L1_1(k+1),'-b');
hold on
plot(k,L1_2_5(k+1),'-r');
hold on
plot(k,L1_3_5(k+1),'-g');
title('Lipschitz constant along iterations, L_0=5');
xlabel('iterations');
ylabel('L_k');
legend('GFISTA (known L_f)','GFISTA cl. backtracking','GFISTA full backtracking');
%%
figure(6)
subplot(2,3,1);
k=0:20;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_20(k+1),'-r');
hold on
plot(k,e1_3_20(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=20');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,2);
k=20:40;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_20(k+1),'-r');
hold on
plot(k,e1_3_20(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=20');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,3);
k=40:60;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_20(k+1),'-r');
hold on
plot(k,e1_3_20(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=20');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,4);
k=60:80;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_20(k+1),'-r');
hold on
plot(k,e1_3_20(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=20');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%----------------------------------------------
subplot(2,3,5);
k=80:iter2;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_20(k+1),'-r');
hold on
plot(k,e1_3_20(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=20');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%------------------------------------------------
subplot(2,3,6);
k=0:iter2;
plot(k,e1_1(k+1),'-b');
hold on
plot(k,e1_2_20(k+1),'-r');
hold on
plot(k,e1_3_20(k+1),'-g');
hold on
plot(k,e1_4(k+1),'-c');
title('Rate of convergence, L_0=20');
xlabel('iterations');
ylabel('F(u^k)-F(u^ * )');
legend('GFISTA (known L_f)','GFISTA cl. backtracking',...
    'GFISTA full backtracking','FISTA (\mu =0)');
%%
figure(7)
k=0:iter2;
plot(k,L1_1(k+1),'-b');
hold on
plot(k,L1_2_20(k+1),'-r');
hold on
plot(k,L1_3_20(k+1),'-g');
title('Lipschitz constant along iterations, L_0=20');
xlabel('iterations');
ylabel('L_k');
legend('GFISTA (known L_f)','GFISTA cl. backtracking','GFISTA full backtracking');
%%
%saving Figurs
if ~exist('output_figures', 'dir')
    mkdir('output_figures')
end
%--------
for i = 1:7
    saveas(i, fullfile('output_figures', sprintf('shape_%d.png', i)))
end