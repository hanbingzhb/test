%% VX公众号：Matlab techniques出品，谨防假冒！
clc;close all
warning off
luoju=55e-2; % 螺距
k=luoju/2/pi; % 螺线方程的系数 r=k theta
L1=341e-2;
D1=L1-27.5e-2*2; % 龙头把手两个孔之间的距离
L2=220e-2;
D2=L2-27.5e-2*2; % 其他凳子把手两个孔之间的距离

%% 先画出部分螺线
theta=16*2*pi:-0.01:0*pi;
r=k*theta;
x=r.*cos(theta);
y=r.*sin(theta);
figure(1)
set(gcf,'Position',[200 200 600 600]);
plot(x,y,'--')
axis equal
grid on
xlabel('x')
ylabel('y')
hold on
% 第一步，先确定300s内，龙头第一个把手的位置演化
mydtheta=@(t,theta)-1./(k*sqrt(1+theta.^2));
theta0=57.032076651015522; % 初始位置时候的角度,我们假设从300秒开始出发，这个角度是一问里面第300秒处的头把手的角度数据！
dt=0.2; % 时间步长
flag=0;
step=0;

N=223; % 龙头+龙身+龙尾总的个数
X=nan*zeros(N+1,3); % 记录每个把手点在一个时间区间内的值
Y=nan*zeros(N+1,3); % 每一行代表每个凳子的前把手孔的位置在各个时间点处的值,因为尾部还有一个孔，所以一共223+1个，所以X和Y对应223+1行
Theta=nan*zeros(N+1,3); % 记录每个孔在时间区间的位置对应的角度theta（这为了求速度的）
Theta(1,3)=theta0; % 这是头把手初始时刻(t=300)已知的角度值

while flag==0 % 如果flag=0，说明一直没有接触，可以继续往前进
    step=step+1;
    X(:,1)=X(:,3);
    Y(:,1)=Y(:,3);
    Theta(:,1)=Theta(:,3); % 开始下一时间步之前，把当前时间区间末尾处的各个把手位移值当作下一个区间开始的数据，所以进行这样的赋值
    tspan=[0,dt/2,dt]; % 求解下一个步长下的位置点
    [tt,theta]=ode45(mydtheta,tspan,Theta(1,1)); % 龙格库塔法求解
    X1=k*theta.*cos(theta);
    Y1=k*theta.*sin(theta); % 下一个时刻，头部把手所在位置坐标
%     for i=1:10:length(theta)
%         title({['t=',num2str(tt(i))],'VX公众号Matlab techniques出品','头部第一个把手中心的轨迹'})
%         plot(X1(i),Y1(i),'r.','MarkerSize',10)
%         drawnow
%     end
%     hwait=waitbar(0,'计算开始...')
    %% 第二步,确定每个时间点下，头部凳子的后面一个孔（也要在螺线上）,以及龙身和龙尾凳子各个孔所在位置(都要在螺线上）
    X(1,:)=X1;
    Y(1,:)=Y1;% 第一行已知了，上面ode45求得的头部第一个把手位置数据
    Theta(1,:)=theta; % 第一行已知了，第一个把手的角度数据，上面求了
    for j=2
        for i=2:N+1 % 在下一个时间点下,求出此时各个把手孔的位置信息
            d=D1*(i<=2)+D2*(i>2); % 分辨下是第一个凳子还是其他凳子，孔之间的距离不一样！
            thetaij=solve_theta(luoju,X(i-1,j),Y(i-1,j),Theta(i-1,j),d); % 子函数求解下一个孔的角度值
            if thetaij>16*2*pi % 如果下一个孔的角度超出了初始位置咋办？不记录！
                Theta(i,j)=nan;
                X(i,j)=nan;
                Y(i,j)=nan;
                break;
            else
                Theta(i,j)=thetaij;
                X(i,j)=k*thetaij*cos(thetaij);
                Y(i,j)=k*thetaij*sin(thetaij);
            end
%             waitbar(((j-1)*N+i)/(length(tt)*N),hwait,'已经完成...')
        end
    end
    hp=plot(X(:,2),Y(:,2),'k-','LineWidth',1.2,'Marker','o','MarkerSize',6,'MarkerFaceColor','r'); % 更新下此时的龙
    title({['t=',num2str(300+step*dt)],'VX公众号Matlab techniques出品','继续盘入，轨迹'})
    drawnow
    % 以上，计算完成了下一个时间点下，各个位置点 下面要判断是否接触！公众号Matlab techniques发布，其他出处皆为抄袭！
    for i=1:N
        x_1=X(i,2);x_2=X(i+1,2);
        y_1=Y(i,2);y_2=Y(i+1,2); % 意思是，对于第一块板，先确定此时两个孔的坐标
        theta_1=Theta(i,2);
        theta_2=Theta(i+1,2);  % 分别是这两个孔的角度
        if isnan(theta_2) % 如果判断到了初始盘入点，剩余的点就不用判断了
            break;
        end
        index1=find((theta_1+2*pi-Theta(:,2))>0); % 从外面一层里面找
        index1=index1(end-2:end); % 找到距离前把手最近的外面一层的三个点的指标
        index2=find(Theta(:,2)-(theta_2+2*pi)>0); 
        if isempty(index2)
            break; % 如果到达了盘入口附近，直接不用判断了，开始下一个时刻位置更新
        else
            index2=index2(1:min(3,length(index2)));     % 找到距离后把手最近的外面一层的三个点的指标 (可能找不到三个点，因为快到初始盘入点时，上面找不到三个最近点了）
        end
        index_i=index1(1):index2(end); % 这是当前全部要考虑的外面一层的把手点的位置指标,从小往大（逆时针排布)
        n=10;m=20; % 板的均匀离散点数量,太密精度高但是运算慢
        for kk=1:length(index_i)-1 % 开始判断第i块板是否与临近的几块板接触
            X2_1=[X(index_i(kk),2);Y(index_i(kk),2)]; % 这是第index_i(kk)把手点所在坐标
            X2_2=[X(index_i(kk+1),2);Y(index_i(kk+1),2)]; % 这是第index_i(kk+1)把手点所在坐标，它两形成第二个长方形，下面就是要判断第i块板是否与这个玩意相交
            panduan=find_if_intersect(L1*(i<=1)+L2*(i>1),[x_1;y_1],[x_2;y_2],L2,X2_1,X2_2,n,m); % 调用子函数，判断是否相交，如果是，返回非空集合panduan，若没有相交，panduan为空集
            if ~isempty(panduan)
                flag=1; % 如果任何一个有相交，设置flag为1,不用再判断了
                break;
            end
        end
        if flag==1  % 如果判断哪一块有接触，直接不用继续下一个时刻了，这个时刻就是最大盘入时刻
            break;
        end
    end
    
    step
    delete(hp)
%     close(hwait)
end
%%
hp=plot(X(:,2),Y(:,2),'k-','LineWidth',1.2,'Marker','o','MarkerSize',6,'MarkerFaceColor','r'); % 更新下此时的龙
title({['t=',num2str(300+step*dt)],'VX公众号Matlab techniques出品','继续盘入，轨迹'})
drawnow



%% 下面来求解速度数据(有各个孔的角度数据Theta，速度数据可以用数值差分方法
V=zeros(size(X)); % 每一行对应一个孔的速度在不同时间节点上的值
V(:,1)=-k*sqrt(1+Theta(:,1).^2).*(Theta(:,2)-Theta(:,1))/dt; % 第一个时间点，用前向差分求导数
V(:,end)=-k*sqrt(1+Theta(:,end).^2).*(Theta(:,end)-Theta(:,end-1))/dt; % 最后一个时间点，用后向差分求导数
V(:,2:end-1)=-k*sqrt(1+Theta(:,2:end-1).^2).*(Theta(:,3:end)-Theta(:,1:end-2))/2/dt; % 中间的每一个时间点，用中心差分来求导数
% 不明白为什么这样求导的，去看我的推文：https://mp.weixin.qq.com/s/rI9B-VHT9i2jA28CN7oFMQ
% 下面验证一下第一个孔位置的速度是不是1，也就是V矩阵的第一行的值
figure
plot(tt,V(1,:),'b-','LineWidth',1.3)
ylim([0 1.1])
xlabel('时间')
ylabel('头把手的速度')
title('验证数值计算得到的速度')
text(100,0.6,'公众号：Matlab techniques')

%% %%%%%%%%%%%%%%%%%%%% 到此，我们已经得到了所有把手中心的位置和速度随时间的数据,下面填表
nn=1/dt;
index=1:nn:length(tt); % 找出特定的时间点,来记录
Dataxy=zeros(2*(N+1),length(index)); % 这些记录所有点的位置信息
Dataxy(1:2:end,:)=round(X(:,index),6); % 这是所有点的x坐标
Dataxy(2:2:end,:)=round(Y(:,index),6); %             y坐标
Datav=round(V(:,index),6); % 这些记录所有点的速度信息
% 写入文件
filename = 'result1_test.xlsx';
sheet = 1;
xlRange = 'B2';
xlswrite(filename,Dataxy,sheet,xlRange)  % 自己搜一下如何保留6位小数（我没管这个细节
sheet = 2;
xlRange = 'B2';
xlswrite(filename,Datav,sheet,xlRange)  % 自己搜一下如何保留6位小数（我没管这个细节

%% 下面再看论文中的格式怎么导出
nn2=60/dt;
index2=1:nn2:length(tt); % 找出特定的时间点,来记录
index_row=[1 2:50:224 224];
Dataxy2=zeros(2*length(index_row),length(index2)); % 这些记录所有点的位置信息
Dataxy2(1:2:end,:)=round(X(index_row,index2),6); % 这是所有点的x坐标
Dataxy2(2:2:end,:)=round(Y(index_row,index2),6); %             y坐标
Datav2=round(V(index_row,index2),6); % 这些记录所有点的速度信息

%%%%%%%%% 以上全部计算完成！%%%%%%%%%%%%
%%%%%%公众号Matlab techniques出品！其他出处皆为抄袭！%%%%%%%%%%%%%%

%% 定义子函数：用来求解螺线上任意一个孔所在点坐标知道的情况下，如何求下一个孔所在位置
function theta=solve_theta(luoju,x1,y1,theta1,d) % 知道等距螺线上任意一点x1,y1,求在这同一条螺线上且与(x1,y1)这点相距为d的点的角度,要求这个角度一定要大于(x1,y1)这点的角度theta1
k=luoju/2/pi;
fun=@(theta)(k*theta.*cos(theta)-x1).^2+(k*theta.*sin(theta)-y1).^2-d^2;% 利用距离
q=0.01;
options = optimoptions('fsolve','Display','off'); % 不提示结果
theta=fsolve(fun,theta1+q,options); % 以前一个孔对应的角度值为基准进行非线性方程求零点，为什么加0.1？很简单，不想找比它还小的角度
while theta<=theta1 || abs(k*theta-k*theta1)>luoju/2 % 如果求解得到的角度比前面的孔对应的角度还小，或者新求出的孔与前一个孔不在一条螺线（在其他螺线）
    q=q+0.1;
    theta=fsolve(fun,theta+q,options); % 重新求一个角度
end % 直到满足：这个点与前面孔在一条螺线，且这个点角度比前面孔的角度大!
end

