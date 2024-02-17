% %Q1.A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
% clear all;
% t=0:2*pi/500:4*pi;
% 
% x=exp(4j*pi*t)+exp(-j*4*pi*t)
% 
% figure(1)
% plot(t,x)
% 
% x=(exp(j*4*pi*t)-exp(-j*4*pi*t))/j;
% 
% figure(2)
% plot(t,x)
% %Q1.B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
% clear all;
% t=0:0.01:8;
% 
% theta=-pi/2;
% z=8*cos(2*pi*t+theta);
% 
% figure(3)
% subplot(4,1,1)
% plot(t,z)
% 
% theta=-pi;
% z=8*cos(2*pi*t+theta);
% 
% figure(3)
% subplot(4,1,2)
% plot(t,z)
% 
% theta=-3*pi/2;
% z=8*cos(2*pi*t+theta);
% 
% figure(3)
% subplot(4,1,3)
% plot(t,z)
% %Q1.C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
% clear all;
% t=0:0.01:8;
% 
% theta=-pi/2;
% w=3*cos(8*pi*t+theta);
% 
% figure(4)
% subplot(4,1,1)
% plot(t,w)
% 
% theta=-pi;
% w=3*cos(8*pi*t+theta);
% 
% figure(4)
% subplot(4,1,2)
% plot(t,w)
% 
% theta=-3*pi/2;
% w=3*cos(8*pi*t+theta);
% 
% figure(4)
% subplot(4,1,3)
% plot(t,w)
% %Q1.D-B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
% clear all;
% t=0:0.01:8;
% 
% T=1;
% tor=0;
% x=8*cos(2*pi*(t-tor/T))
% 
% figure(5)
% subplot(4,1,1)
% plot(t,x)
% grid on 
% 
% tor=0.25;
% x=8*cos(2*pi*(t-tor/T))
% subplot(4,1,2)
% plot(t,x)
% grid on
% 
% tor=0.75;
% x=8*cos(2*pi*(t-tor/T))
% subplot(4,1,3)
% plot(t,x)
% grid on
% % Q1.D-C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
% clear all;
% t=0:0.01:8;
% T=0.25;
% tor=0;
% x=3*cos(2*pi*(t-tor/T));
% 
% figure(6)
% subplot(4,1,1)
% plot(t,x)
% grid on
% 
% tor=0.0625;
% x=3*cos(2*pi*(t-tor/T));
% 
% figure(6)
% subplot(4,1,2)
% plot(t,x)
% 
% grid on
% tor=0.125;
% x=3*cos(2*pi*(t-tor/T));
% 
% figure(6)
% subplot(4,1,3)
% plot(t,x)
% grid on
% % % Q2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Example code
% clear, clc
% n=100; %define the number of terms of the F.S>
% %define the time index. The range is from[0:1],but with
% %step size of 0.002 for fine resolution
% t = 0:0.002:2;
% T=2;
% for i = 1:n
%     x(i,:) =(4/pi)*(1/((2*i)-1)*sin(2*pi*((2*i)-1)*t/T));
% end
% %plot the first term, the second term
% %the third term, and the Sum of teh first three terms
% %approch 1:
% figure(7)
% plot(t,x(1,:)),hold on
% plot(t,x(2,:),'g')
% plot(t,x(3,:),'m')
% plot(t,x(1,:)+x(2,:)+x(3,:),'r')
% hold off
% xlabel('t')
% legend('x(1)','x(2)','x(3)','x(1)+x(2)+x(3)')
% %approch2:
% figure(7)
% subplot(5,1,1)
% plot(t,x(1,:)) %plot the first term
% subplot(5,1,2)
% plot(t,x(2,:),'g')%plot the second term
% subplot(5,1,3)
% plot(t,x(3,:),'m')%plot the third term
% subplot(5,1,4)
% plot(t,x(4,:),'m')%plot the third term
% subplot(5,1,5)
% plot(t,sum(x(1:n,:),1),'r')
% %Q2.A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear, clc
% N = 100; % define the number of terms of the F.S.
% t = 0:0.002:2; % define the time index with fine resolution
% T = 2;
% omega_0 = 2*pi/T;%!!!!!!!this is omega value!!!!!!!!!!!!it valiue is PI!!!!!!
% x_N = zeros(size(t));
% for n = 1:N
%     x_N = x_N + (4/pi)*(1/(2*n - 1))*sin((2*n - 1)*omega_0*t);
% end
% figure(9)
% subplot(5,1,1)
% plot(t, (4/pi)*(1/(2*1 - 1))*sin((2*1 - 1)*omega_0*t), 'b')
% title('x(1)')
% subplot(5,1,2)
% plot(t, (4/pi)*(1/(2*2 - 1))*sin((2*2 - 1)*omega_0*t), 'g')
% title('x(2)')
% subplot(5,1,3)
% plot(t, (4/pi)*(1/(2*3 - 1))*sin((2*3 - 1)*omega_0*t), 'm')
% title('x(3)')
% subplot(5,1,4)
% plot(t, x_N, 'r')
% title('x(1) + x(2) + x(3)')
% subplot(5,1,5)
% plot(t, (4/pi)*(1/(2*4 - 1))*sin((2*4 - 1)*omega_0*t), 'k')
% title('x(4)')
% xlabel('t')
% %Q2.B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear, clc
% n = 500; % Set the number of terms to 50
% t = 0:0.002:2; % Time vector
% T = 2;
% 
% % Initialize x to store the terms
% x = zeros(n, length(t));
% 
% for i = 1:n
%     x(i,:) = (4/pi)*(1/((2*i)-1)*sin(2*pi*((2*i)-1)*t/T));
% end
% 
% % Plot individual terms and their sum
% figure(10)
% subplot(6,1,1:4)
% for i = 1:4
%     plot(t, x(i,:), 'DisplayName', ['x(', num2str(i), ')']), hold on
% end
% hold off
% xlabel('t')
% legend('show')
% 
% subplot(6,1,5:6)
% plot(t, sum(x(1:n,:), 1), 'r', 'DisplayName', 'Approximation (N=500)')
% xlabel('t')
% legend('show')
% %Q2.C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Example code
% clear, clc
% n=40; %define the number of terms of the F.S>
% %define the time index. The range is from[0:1],but with
% %step size of 0.002 for fine resolution
% t = 0:0.002:2;
% T=2;
% for i = 1:n
%     x(i,:) =(4/pi)*(1/((2*i)-1)*sin(2*pi*((2*i)-1)*t/T));
% end
% %plot the first term, the second term
% %the third term, and the Sum of teh first three terms
% %approch:
% figure(12)
% subplot(5,1,1)
% plot(t,x(1,:)) %plot the first term
% subplot(5,1,2)
% plot(t,x(2,:),'g')%plot the second term
% subplot(5,1,3)
% plot(t,x(3,:),'m')%plot the third term
% subplot(5,1,4)
% plot(t,x(4,:),'m')%plot the third term
% subplot(5,1,5)
% plot(t,sum(x(1:n,:),1),'r')
% hold on 
% unitstep1 = t>=0; %equavalent to if t >= 0 unitstep1 = True else unitstep1 =False end
% unitstep2 = t>=1;
% figure (13)
% plot(t,unitstep1-2*unitstep2)
% hold off
% %Q2.D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear, clc
% 
% % Define parameters
% t = 0:0.002:2; % Define the time index with fine resolution
% T = 2;
% omega_0 = 2*pi/T;
% 
% % Generate the square wave function x(t)
% x_t = square(2*pi*t/T);
% 
% % Initialize variables for different approximations
% N_values = [10, 100, 500];
% num_N_values = numel(N_values);
% 
% % Plot subplots
% figure(8)
% 
% for i = 1:num_N_values
%     N = N_values(i);
% 
%     % Calculate Fourier series approximation
%     x_N = zeros(size(t));
% 
%     for n = 1:N
%         x_N = x_N + (4/pi)*(1/(2*n - 1))*sin((2*n - 1)*omega_0*t);
%     end
% 
%     % Calculate the approximation error e_a(t)
%     e_a_t = x_t - x_N;
% 
%     % Plot results in subplots
%     subplot(num_N_values, 2, 2*i-1)
%     plot(t, x_N, 'r', 'DisplayName', ['\hat x_{' num2str(N) '}(t)'])
%     xlabel('t')
%     ylabel(['\hat x_{' num2str(N) '}(t)'])
%     title(['Approximation (N=' num2str(N) ')'])
%     legend('show')
% 
%     subplot(num_N_values, 2, 2*i)
%     plot(t, e_a_t, 'g', 'DisplayName', ['Error (N=' num2str(N) ')'])
%     xlabel('t')
%     ylabel(['e_a_{' num2str(N) '}(t)'])
%     title(['Approximation Error (N=' num2str(N) ')'])
%     legend('show')
% end
% %Q2.E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear, clc
% 
% % Define parameters
% t = 0:0.002:2; % Time vector
% T = 2;
% num_terms = 20; % Maximum number of terms to consider
% 
% mse_values = zeros(1, num_terms);
% for n = 1:num_terms
%     % Compute Fourier Series terms
%     x = zeros(n, length(t));
%     for i = 1:n
%         x(i,:) = (4/pi)*(1/((2*i)-1)*sin(2*pi*((2*i)-1)*t/T));
%     end
% 
%     % Compute the sum of the first n terms
%     syw = sum(x(1:n,:), 1);
% 
%     % Generate square wave
%     unitstep1 = t >= 0;
%     unitstep2 = t >= 1;
%     sqw = unitstep1 - 2*unitstep2;
% 
%     % Compute Mean Square Error (MSE)
%     er = mean((sqw - syw).^2);
%     mse_values(n) = er;
% end
% 
% % Plot the Mean Square Error (MSE) with regard to the number of terms
% figure(15)
% plot(1:num_terms, mse_values, 'bo-', 'LineWidth', 2)
% xlabel('Number of Terms')
% ylabel('Mean Square Error (MSE)')
% title('Convergence of Fourier Series Approximation')
% grid on









