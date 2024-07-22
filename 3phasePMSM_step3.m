clear all; clc; close all;
%입력 자기장은 1, 전기각에 대한 한 주기 2, 사인파의 형태로 입력을 해줘야함(입력전류가 사인파로 정의 되었기 때문)

%%
%변수선언

D = 0.046; %회전자 지름
L = 0.07; %적층 높이
N = 48*4; %턴 수
RPM = 2000; %RPM
p = 8; %자석의 개수
max_harmonics = 13; %파악하는 고조파의 개수

theta_w = pi*0.6; %short pitch angle

W_factor = 0; 
%distribution factor(q, theta_d)를 고려하면 1, 고려하지 않으면 0 
% W_factor = 0 이면 q, theta_d 의 값은 사용되지 않음. 
q = 1; %number of distributed coils
theta_d = 0; %distribution angle

%%%%%%%%%%%%%%%%%%%%수정금지%%%%%%%%%%%%%%%%%%%%%%%
N_fft = 201; %fft 샘플링 갯수(=전기각 해상도) (N_fft > 2*max_harmonics 의 조건을 만족해야함, 3의 배수여야 함)
w_m = RPM/60*2*pi; %기계각속도
w_e = w_m*p/2; %전기각속도
f_m = w_m/(2*pi); %기계주파수
f_e = w_e/(2*pi); %전기주파수
n_theta_e = N_fft; %전기각 해상도
theta_e = linspace(0, 2*pi,n_theta_e); %전기각
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_txt = readmatrix('B_data.txt'); %0.785*sqrt(1)*square(theta_e)%전기각에 대한 자기장
B_theta_e = B_txt;
n_sample = length(B_theta_e); %자기장 데이터의 길이
B_theta_e = resample(B_theta_e, N_fft, n_sample); %자기장의 데이터 N_fft(3의 배수)로 resample
%3상으로 자기장 변환
%B2_theta_e = zeros(1, N_fft);
%B2_theta_e(1:N_fft/3) = B_theta_e(N_fft/3*2+1:end);
%B2_theta_e(N_fft/3+1:end) = B_theta_e(1:N_fft/3*2);
%B3_theta_e = zeros(1, N_fft);
%B3_theta_e(1:N_fft/3) = B2_theta_e(N_fft/3*2+1:end);
%B3_theta_e(N_fft/3+1:end) = B2_theta_e(1:N_fft/3*2);

%전기각에 대한 입력전류
I_amp = 14*sqrt(1); %입력전류의 진폭
I1_theta_e = I_amp*sin(theta_e);
I2_theta_e = I_amp*sin(theta_e-(pi*2/3));
I3_theta_e = I_amp*sin(theta_e-(pi*4/3));

%%
%fft를 이용하여 B_t의 푸리에 급수 계수 구하기

f_fft = f_e*(0:N_fft-1); %이산주파수(샘플링 주파수가 f_e*N_fft 이다. > 주파수 도메인 1칸 주파수는 f_e 임)

B_theta_e_fft = fft(B_theta_e); %B_theta_e에 대한 fft

cutoff_positive = ceil(N_fft/2); %양수에 해당하는 N_fft 범위
cutoff_max_harmonics = max_harmonics+1; %설정한 고조파 수에 해당하는 N_fft 범위

B_theta_e_DTFS = B_theta_e_fft(1:cutoff_positive)/N_fft; %B_t의 DTF를 DTFS로 변환
B_theta_e_DTFS_coefficient = B_theta_e_DTFS(1:cutoff_max_harmonics); %설정한 고조파의 수까지 자름
B_theta_e_DTFS_coefficient(2:cutoff_max_harmonics)=2*B_theta_e_DTFS_coefficient(2:cutoff_max_harmonics); %푸리에급수 계수로 변환
f_DTFS_coefficient = f_fft(1:cutoff_max_harmonics); %B_t_DTFS_coefficient의 그래프 plot을 위한 이산주파수 벡터 생성

B_g = zeros(1,max_harmonics); %자기장 푸리에 급수의 계수

for i = 1:max_harmonics %고조파 샘플링
    B_g(i) = max(abs(B_theta_e_DTFS_coefficient(i+1)));
end


%%
%역기전력의 계산

K_e = D*L*N*B_g; %Ke의 계산
E = K_e*w_m; %고조파의 계수
K_w = zeros(1, max_harmonics);

for i = 1:max_harmonics %winding factor(단절계수, 분포계수)고려
    if W_factor == 1 %winding factor
        K_w(i) = cos(i/2*theta_w)*( sin(i*q*theta_d/2)/q/sin(i*theta_d/2) );
    else
        K_w(i) = cos(i/2*theta_w);
    end
    E(i) = E(i)*K_w(i);
end

e_theta_e = zeros(max_harmonics, n_theta_e); %전기각에 대한 역기전력 고조파
e_theta_e_sum = zeros(1, n_theta_e); %sum

for i = 1:max_harmonics %역기전력 고조파
    e_theta_e(i,:) = E(i)*sin(theta_e*i);
end

for i = 1:max_harmonics %역기전력 고조파를 모두 sum
    e_theta_e_sum = e_theta_e_sum + e_theta_e(i,:);
end

%%
%출력 계산

P_theta_e = zeros(max_harmonics, n_theta_e); %전기각에 대한 파워 고조파
P_theta_e_sum = zeros(1, n_theta_e); %sum

for i = 1:max_harmonics %3상의 power 고조파
    P_theta_e(i,:) = E(i)*sin(theta_e*i).*I1_theta_e + E(i)*sin((theta_e-(pi*2/3))*i).*I2_theta_e + E(i)*sin((theta_e-(pi*4/3))*i).*I3_theta_e;
end

for i = 1:max_harmonics %power의 고조파를 모두 sum
    P_theta_e_sum = P_theta_e_sum + P_theta_e(i,:);
end

%%
%토크계산

T_theta_e = P_theta_e/w_m; %전기각에 대한 토크 고조파
T_theta_e_sum = zeros(1, n_theta_e); %sum

for i = 1:max_harmonics %토크의 고조파를 모두 sum
    T_theta_e_sum = T_theta_e_sum + T_theta_e(i,:);
end


%%
%그래프

figure(1);
tiledlayout('vertical');
nexttile;
plot(theta_e, B_theta_e);
title('공극에서 전기각 한 주기의 자기장')
xlabel('전기각(rad)')
ylabel('자기장(T)')

nexttile;
plot(f_DTFS_coefficient, abs(B_theta_e_DTFS_coefficient), 'o');
title('자기장 고조파의 계수')
xlabel('frequency(Hz)')
ylabel('자기장 고조파의 크기(T)')

figure(2);
plot(theta_e, e_theta_e_sum);
title('전기각에 대한 역기전력')
xlabel('전기각(rad)')
ylabel('역기전력(V)')

figure(3);
plot(theta_e, P_theta_e_sum);
title('전기각에 대한 출력')
xlabel('전기각(rad')
ylabel('출력(VI=W)')

figure(4)
tiledlayout('vertical');
nexttile;
plot(theta_e, T_theta_e_sum);
title('전기각에 대한 토크(sum)')
xlabel('전기각(rad)')
ylabel('토크(N*m)')

nexttile;
plot(theta_e, T_theta_e(1,:)+T_theta_e(3,:)+T_theta_e(9,:));
title('1, 3, 9고조파 토크')
xlabel('전기각(rad)')
ylabel('토크(N*m)')

nexttile;
plot(theta_e, T_theta_e(5,:)+T_theta_e(6,:)+T_theta_e(7,:));
title('6고조파 토크 리플')
xlabel('전기각(rad)')
ylabel('토크(N*m)')

nexttile;
plot(theta_e, T_theta_e(11,:)+T_theta_e(12,:)+T_theta_e(13,:));
title('12고조파 토크 리플')
xlabel('전기각(rad)')
ylabel('토크(N*m)')
