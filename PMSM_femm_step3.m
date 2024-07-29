clear all; clc; close all;
addpath("C:/femm42/mfiles/");

%%
%<변수선언>

L = 70; %적층높이
D = 46; %직경
D_2 = 42; %중심 ferrite의 직경
p = 8; %자석의 개수
a_p = 0.89; %극호비
w_m = D_2*pi/p*a_p; %자석의 너비(호의 길이)
h_m = 1; %자석의 높이
g =(D-D_2)/2-h_m; %공극 길이(stator 반지름 - 중심 ferrite 반지름)
arc_m = 360/p*1.5; %자석의 아크(deg)
num_slot = 12; %slot의 개수
theta_slot = pi/10; %하나의 slot이 기계각에서 차지하는 각(rad)
b_so = 0.5; %slot 오프닝 폭(mm)
background = 100; %도화지의 가로, 세로 길이(중심은 원점)
N = 200; %샘플링 개수
B_txt = 0; %B_data.txt를 만들면 1, 만들지 않으면 0

%%%%%%%%%%%%%%%%%%%%카터계수%%%%%%%%%%%%%%%%%%%%%%
Cater = 1; %1이면 카터계수가 고려된 모델을 그림, 아니면 실제의 형상을 그림

if Cater == 1
    tau = pi*(D+2*g)/num_slot;
    beta = 4/pi*(b_so/(2*(h_m+g))*atan(b_so/(2*(h_m+g)))-log(sqrt(1+(b_so/(2*(h_m+g)))^2)));
    k_c = tau/(tau-beta*(h_m+g)); %카터계수
    g_e = g+(k_c-1)*(h_m+g); %카터계수가 고려된 에어갭
    delta_D = g_e-g; %D의 변화

    D = D+delta_D;
    b_so = 0; %slot 오프닝 폭 초기화
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_s = D+30; % stator 직경
D_c1 = D+2; % slot 안쪽 점 직경
D_c2 = D+22; % slot  바깥 점 직경


%%
%<기초 설정>

openfemm; %FEMM실행
newdocument(0); %새 Document 창 실행(정자계 해석)

%파일 저장
path ='C:/femm42/'; %Path 설정
name_fem ='spm_motor.fem'; %파일 이름 설정
mi_saveas([path, name_fem]); %파일 저장

%재료 불러오기
Hc=900000;
mi_addmaterial('N4X', 1.05, 1.05, Hc, 0, 0.56); %mi_addmaterial("materialname", mu x, mu y, H c, J, Cduct, Lam d, Phi hmax, lam fill, LamType, Phi hx, Phi hy,NStrands,WireD)
mi_getmaterial('Air'); 
mi_getmaterial('M-19 Steel');
mi_modifymaterial('M-19 Steel', 0, '35H270'); % mi_modifymaterial("BlockName",propnum,value)
mi_addmaterial('a+', 1, 1, 0, 0, 56)
mi_addmaterial('a-', 1, 1, 0, 0, 56)
mi_addmaterial('b+', 1, 1, 0, 0, 56)
mi_addmaterial('b-', 1, 1, 0, 0, 56)
mi_addmaterial('c+', 1, 1, 0, 0, 56)
mi_addmaterial('c-', 1, 1, 0, 0, 56)
mi_addboundprop('A_0', 0, 0, 0, 0);	%도화지의 태두리 기본 설정 %mi_addboundprop("propname", A0, A1, A2, Phi, Mu, Sig, c0, c1, BdryFormat, ia, oa)

%도화지 설정
depth=L; % 화면방향 깊이
mi_probdef(0,'millimeters','planar', 1e-009, depth, 30) % miilimeters = 단위 사용, planar = 평면 사용 %mi_probdef(frequency,units,type,precision,(depth),(minangle),(acsolver))
drawrectangle(-background/2, -background/2, background/2, background/2);
selectrectangle(-background/2, -background/2, background/2, background/2); 
mi_setsegmentprop('A_0', 0, 1, 0, 100) % 바운더리 컨디션 설정 % mi_setsegmentprop("propname", elementsize, automesh, hide, group)
mi_clearselected()

mi_zoomnatural() %전체화면으로 맞춤

%%
%<밑그림 그리기>

drawcircle(0, 0, D/2); %직경 그리기
drawcircle(0, 0, D_2/2); %중심 ferrite
drawcircle(0, 0, D_s/2) %stator 그리기

%%
%<자석 한 개 그리기>

theta_m = 2*pi/p; %자석이 서로 떨어진 각

%자석의 점
m_dot = zeros(2, 4);
m_dot(:, 1) = (D_2/2)*[cos(0); sin(0)];
m_dot(:, 2) = (D_2/2+h_m)*[cos(w_m/(D_2/2)/2); sin(w_m/(D_2/2)/2)];
m_dot(:, 3) = (D_2/2)*[cos(w_m/(D_2/2)); sin(w_m/(D_2/2))];
m_dot(:, 4) = (D_2/2+h_m/2)*[cos(w_m/(D_2/2)/2); sin(w_m/(D_2/2)/2)]; %label dot
mi_addnode(m_dot(1, 1), m_dot(2, 1));
mi_addnode(m_dot(1, 2), m_dot(2, 2));
mi_addnode(m_dot(1, 3), m_dot(2, 3));

%자석 그리기
mi_addarc(m_dot(1, 1), m_dot(2, 1), m_dot(1, 2), m_dot(2, 2), arc_m/2, 5);
mi_addarc(m_dot(1, 2), m_dot(2, 2), m_dot(1, 3), m_dot(2, 3), arc_m/2, 5);

%%
%<slot 한 개 그리기>

theta_s = 2*pi/num_slot; %slot이 서로 떨어진 각

%slot의 점
s_dot = zeros(2, 14);
s_dot(:, 1) = (D_c2/2)*[cos(0); sin(0)];
s_dot(:, 2) = (D_c2/2)*[cos(theta_slot/2); sin(theta_slot/2)];
s_dot(:, 3) = (D_c2/2)*[cos(theta_slot); sin(theta_slot)];
s_dot(:, 4) = (D_c1/2)*[cos(theta_slot/10*9); sin(theta_slot/10*9)];
s_dot(:, 5) = (D/2)*[cos(theta_slot/2+b_so/2/(D/2)); sin(theta_slot/2+b_so/2/(D/2))];
s_dot(:, 6) = (D/2)*[cos(theta_slot/2-b_so/2/(D/2)); sin(theta_slot/2-b_so/2/(D/2))];
s_dot(:, 7) = (D_c1/2)*[cos(theta_slot/10); sin(theta_slot/10)];
s_dot(:, 8) = (D/2)*[cos(theta_slot/2); sin(theta_slot/2)];
s_dot(:, 9) = (D_c1+D_c2)/4*[cos(theta_slot/4*3); sin(theta_slot/4*3)]; %label dot a+
s_dot(:, 10) = (D_c1+D_c2)/4*[cos(theta_s+theta_slot/4); sin(theta_s+theta_slot/4)]; %label dot a-
s_dot(:, 11) = (D_c1+D_c2)/4*[cos(theta_s+theta_slot/4*3); sin(theta_s+theta_slot/4*3)]; %label dot b+
s_dot(:, 12) = (D_c1+D_c2)/4*[cos(2*theta_s+theta_slot/4); sin(2*theta_s+theta_slot/4)]; %label dot b-
s_dot(:, 13) = (D_c1+D_c2)/4*[cos(2*theta_s+theta_slot/4*3); sin(2*theta_s+theta_slot/4*3)]; %label dot c+
s_dot(:, 14) = (D_c1+D_c2)/4*[cos(3*theta_s+theta_slot/4); sin(3*theta_s+theta_slot/4)]; %label dot c-

mi_addnode(s_dot(1, 1), s_dot(2, 1));
mi_addnode(s_dot(1, 2), s_dot(2, 2));
mi_addnode(s_dot(1, 3), s_dot(2, 3));
mi_addnode(s_dot(1, 4), s_dot(2, 4));
mi_addnode(s_dot(1, 5), s_dot(2, 5));
mi_addnode(s_dot(1, 6), s_dot(2, 6));
mi_addnode(s_dot(1, 7), s_dot(2, 7));
mi_addnode(s_dot(1, 8), s_dot(2, 8));

%slot 그리기
mi_addarc(s_dot(1, 1), s_dot(2, 1), s_dot(1, 2), s_dot(2, 2) ,theta_slot/2/pi*180 ,5);
mi_addarc(s_dot(1, 2), s_dot(2, 2), s_dot(1, 3), s_dot(2, 3) ,theta_slot/2/pi*180 ,5);
mi_addsegment(s_dot(1, 3), s_dot(2, 3), s_dot(1, 4), s_dot(2, 4));
mi_addsegment(s_dot(1, 4), s_dot(2, 4), s_dot(1, 5), s_dot(2, 5));
mi_addsegment(s_dot(1, 6), s_dot(2, 6), s_dot(1, 7), s_dot(2, 7));
mi_addsegment(s_dot(1, 7), s_dot(2, 7), s_dot(1, 1), s_dot(2, 1));
mi_addsegment(s_dot(1, 2), s_dot(2, 2), s_dot(1, 8), s_dot(2, 8));

%%
%<회전하여 모두 그리기>

transform_m = [cos(theta_m) -sin(theta_m);
               sin(theta_m) cos(theta_m)]; %자석 회전행렬
transform_s = [cos(theta_s) -sin(theta_s);
               sin(theta_s) cos(theta_s)]; %slot 회전변환행렬

%자석
m_dot_bin = m_dot;
for i = 1:p-1
    m_dot_bin = transform_m*m_dot_bin;
    mi_addnode(m_dot_bin(1, 1), m_dot_bin(2, 1));
    mi_addnode(m_dot_bin(1, 2), m_dot_bin(2, 2));
    mi_addnode(m_dot_bin(1, 3), m_dot_bin(2, 3));
    mi_addarc(m_dot_bin(1, 1), m_dot_bin(2, 1), m_dot_bin(1, 2), m_dot_bin(2, 2), arc_m/2, 5);
    mi_addarc(m_dot_bin(1, 2), m_dot_bin(2, 2), m_dot_bin(1, 3), m_dot_bin(2, 3), arc_m/2, 5);
end

%slot
s_dot_bin = s_dot;
for i = 1:num_slot-1
    s_dot_bin = transform_s*s_dot_bin;
    mi_addnode(s_dot_bin(1, 1), s_dot_bin(2, 1));
    mi_addnode(s_dot_bin(1, 2), s_dot_bin(2, 2));
    mi_addnode(s_dot_bin(1, 3), s_dot_bin(2, 3));
    mi_addnode(s_dot_bin(1, 4), s_dot_bin(2, 4));
    mi_addnode(s_dot_bin(1, 5), s_dot_bin(2, 5));
    mi_addnode(s_dot_bin(1, 6), s_dot_bin(2, 6));
    mi_addnode(s_dot_bin(1, 7), s_dot_bin(2, 7));
    mi_addnode(s_dot_bin(1, 8), s_dot_bin(2, 8));
    mi_addarc(s_dot_bin(1, 1), s_dot_bin(2, 1), s_dot_bin(1, 2), s_dot_bin(2, 2) ,theta_slot/2/pi*180 ,5);
    mi_addarc(s_dot_bin(1, 2), s_dot_bin(2, 2), s_dot_bin(1, 3), s_dot_bin(2, 3) ,theta_slot/2/pi*180 ,5);
    mi_addsegment(s_dot_bin(1, 3), s_dot_bin(2, 3), s_dot_bin(1, 4), s_dot_bin(2, 4));
    mi_addsegment(s_dot_bin(1, 4), s_dot_bin(2, 4), s_dot_bin(1, 5), s_dot_bin(2, 5));
    mi_addsegment(s_dot_bin(1, 6), s_dot_bin(2, 6), s_dot_bin(1, 7), s_dot_bin(2, 7));
    mi_addsegment(s_dot_bin(1, 7), s_dot_bin(2, 7), s_dot_bin(1, 1), s_dot_bin(2, 1));
    mi_addsegment(s_dot_bin(1, 2), s_dot_bin(2, 2), s_dot_bin(1, 8), s_dot_bin(2, 8));
end

%%
%<물성치 넣기>

%자석
name_m = 'N4X';
angle_bin = w_m/(D_2/2)/2/pi*180; %극의 방향
m_dot_bin = m_dot;

for i = 1:p
    if rem(i, 2) == 1
        add_label(m_dot_bin(1, 4), m_dot_bin(2, 4), name_m, angle_bin, 1);
    else
        add_label(m_dot_bin(1, 4), m_dot_bin(2, 4), name_m, angle_bin+180, 2);
    end
    m_dot_bin = transform_m*m_dot_bin;
    angle_bin = angle_bin + theta_m/pi*180;
end

%slot
s_dot_bin = s_dot;

transform_s_label = [cos(theta_s*3) -sin(theta_s*3); %slot label 회전행렬
                     sin(theta_s*3) cos(theta_s*3)];

for i = 1:num_slot/3
    add_label(s_dot_bin(1, 9), s_dot_bin(2, 9), 'a+', 0, 3);
    add_label(s_dot_bin(1, 10), s_dot_bin(2, 10), 'a-', 0, 4);
    add_label(s_dot_bin(1, 11), s_dot_bin(2, 11), 'b+', 0, 5);
    add_label(s_dot_bin(1, 12), s_dot_bin(2, 12), 'b-', 0, 6);
    add_label(s_dot_bin(1, 13), s_dot_bin(2, 13), 'c+', 0, 7);
    add_label(s_dot_bin(1, 14), s_dot_bin(2, 14), 'c-', 0, 8);
    s_dot_bin = transform_s_label*s_dot_bin;
end

%ferrite
name_f = '35H270';
add_label(0, 0, name_f, 0, 9);
add_label((D_s+D_c2)/4, 0, name_f, 0, 9);

%air
name_a = 'Air';
add_label((D_2+D)/4, 0, name_a, 0, 10);
add_label(D_s/2+5, 0, name_a, 0, 10);

%%
% 자기장 분석

%컨투어 시작점과 끝점
start_point = D_2*pi/8*(1-a_p)/2; %반 주기에서 자석이 없는 부위의 절반 길이
start_point_theta = -start_point/(D_2/2); %각으로 환산
end_point_theta = start_point_theta + 2*pi*(2/p);
dot_start = (D/2+D_2/2+h_m)/2*[cos(start_point_theta); sin(start_point_theta)];
dot_end = (D/2+D_2/2+h_m)/2*[cos(end_point_theta); sin(end_point_theta)];

%샘플포인트
theta_N = (0:N-1)/N*2*pi*(2/p); %샘플링 포인트의 각
dot_sample = zeros(2, N);
for i = 1:N
    dot_sample(1, i) = (D/2+D_2/2+h_m)/2*cos(start_point_theta+theta_N(i));
    dot_sample(2, i) = (D/2+D_2/2+h_m)/2*sin(start_point_theta+theta_N(i));
    %mi_addnode(dot_sample(1, i), dot_sample(2, i));
end

%%
% FEMM 해석 수행
mi_analyze(1); % 해석 수행
mi_loadsolution(); % 해석 결과 불러오기

B_xyn = zeros(3, N); %row1은 Bx, row2는 By, row3는 B.n
for i = 1:N
    pointvalue = mo_getpointvalues(dot_sample(1, i), dot_sample(2, i));
    B_xyn(1, i) = pointvalue(2);
    B_xyn(2, i) = pointvalue(3);
    theta_N_nvector = [cos(start_point_theta+theta_N(i)) sin(start_point_theta+theta_N(i))];
    B_xyn(3, i) = theta_N_nvector*B_xyn(1:2, i);
end

figure(1);%그래프 plot
plot(theta_N, B_xyn(3,:), "Color", 'r',  'LineWidth', 3);
xlabel('기계각(rad)');
ylabel('자기장의 법선성분(T)');
if B_txt == 1
    writematrix(B_xyn(3, :), 'B_data.txt'); %txt파일로 추출하기
end

%slot의 넓이 측정
mo_selectblock(s_dot(1, 9), s_dot(2, 9));
A_s = mo_blockintegral(5)*2*num_slot*1000*1000; %단위는 mm^2

%%%%%%%%%%%%%%%컨투어로 femm에서 그래프 확인하기%%%%%%%%%%%%%%%%%%%%%
%mo_addcontour(x_start, y_start);
%mo_addcontour(x_end, y_end);
%mo_bendcontour((2*pi*(2/num_pole))/pi*180, 10);
%mo_makeplot(2, N, 'Bn.txt', 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%함수

function drawrectangle(x1,y1,x2,y2) %(x1,y1)과 (x2,y2) 점을 있는 사각형 그리기 함수
    mi_addnode(x1,y1);	mi_addnode(x1,y2);	mi_addnode(x2,y1);	mi_addnode(x2,y2);
    mi_addsegment(x1,y1,x1,y2);	mi_addsegment(x1,y2,x2,y2);	mi_addsegment(x2,y2,x2,y1);	mi_addsegment(x2,y1,x1,y1); % mi_addsegment(x1,y1,x2,y2)
end

function drawcircle(x1,y1,r) %(x1, y1) 점을 기준으로 r 반지름의 원 그리기 함수
    mi_addnode(x1+r, y1);
    mi_addnode(x1-r, y1);
    mi_addarc(x1+r, y1,x1-r, y1, 180, 5); %mi_addarc(x1,y1,x2,y2,angle,maxseg)
    mi_addarc(x1-r, y1, x1+r, y1, 180, 5);
end

function add_label(x1,y1,name,angle,group) %공간에 물성치 넣어주는 함수
    mi_clearselected();
    mi_addblocklabel(x1, y1); % mi_addblocklabel(x,y): Add a new block label at (x,y)
    mi_selectlabel(x1, y1); %Select the label closet to (x,y). Returns the coordinates of the selected label
    mi_setblockprop(name, 1, 0, 'None', angle, group, 0); %mi_setblockprop("blockname", automesh, meshsize, "incircuit", magdirection, group, turns)
end

function selectrectangle(x1,y1,x2,y2)
    mi_selectsegment((x1+x2)/2, y1);	mi_selectsegment((x1+x2)/2, y2);
    mi_selectsegment(x1, (y1+y2)/2);	mi_selectsegment(x2, (y1+y2)/2);
end


