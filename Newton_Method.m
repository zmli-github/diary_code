%定义函数
%{
function Fp = F(p)
Fp = P_0*p^3-3*P_0*(norm(m-p*A,2))^2*p/U_tip^2-6*P_0*A'*(m-p*A)*p^2/U_tip-d_0*p_a*s*B*(norm(m-p*A,2))^3
-'Par_Phi_t'*p^3-eta_A*'laplace_Phi'*p^3+1/tau*(p-p_k)*p^3-1.5*d_0*p_a*s*B*norm(m-p*A,2)*A'*(m-p*A)*p ;%（13）式的表达式
end
function DFp = DF(p)
DFp = 3*P_0*p^2-3*P_0*(-2*A'*(m-p*A)*p+(norm(m-p*A,2))^2)-6*P_0(-A'*A*p^2+2*A*(m-p*A)*p)/U_tip
-d_0*p_a*s*B((-2*A'*(m-p*A))*norm(m-p*A,2)+(norm(m-p*A,2))^2*(-A'*(m-p*A)/norm(m-p*A,2)))-3*'Par_Phi_t'*p^2
-3*eta_A*'laplace_Phi'*p^2+(4*p^3-3*p_k*p^2)/tau-1.5*d_0*p_a*s*B(norm(m-p*A,2)*A'*(m-p*A)
+p*(-A'(m-p*A)*A'*(m-p*A)/norm(m-p*A,2)-norm(m-p*A)*A'*A));%（13）式的导函数表达式
end   
function laplace_p = Laplace_p(p)
laplace_p = 
%   函数体展示\Phi的迭代式
  

end
function par_p_t = PAR_p_t(p)
par_p_t = 1000*(p_k1-p_k)
%   函数体展示\Phi的迭代式
  

end
%}

%主程序
clc;
clear;
p=[1,0.999];m=[[0.01,0.011];[0.01,0.011]];%迭代初始值,'关于这个m1,m2值的问题需要商议'大概率需要调试
Phi=[0.1,0.11];
P_0=1e-7;p_a=1.225;d_0=0.3;s=0.05;B=0.01;U_tip=3.16e-4;
eta_A=1e-2;%设定平均风速的方差为：1e-2,1e-2.5,1e-5
A=input('请输入平均风速（2x1列向量形式）：');
dt=0.0001;tau=0.1;sigma=0.1;

format long
I=input('请输入想要迭代的总次数（整数）');  
i=1;
Vxy=[0.007;0.007]; Dxy=[]; X=[0.1;0.1];
grad_p=[0;0];grad_Phi=[0;0];grad_m=[0;0];
laplace_p=[];laplace_Phi=[];neg_laplace_Phi=[];
while i<I%迭代次数
    
    Vxy(:,i)=m(:,i)./p(i)-A;
    Dxy(:,i)=Vxy(:,i)*dt;
    X(:,i+1)=X(:,i)+Dxy(:,i);
    par_p_t=(p(i+1)-p(i))/dt;
    par_Phi_t=(Phi(i+1)-Phi(i))/dt;
    
    %grad_p(i+1)=[([p(i),X(1,i+1),X(2,i)]-[p(i),X(:,i)])/dt,([p(i),X(1,i),X(2,i+1)]-[p(i),X(:,i)])/dt];
    %grad_Phi(i+1)=[([Phi(i),X(1,i+1),X(2,i)]-[Phi(i),X(:,i)])/dt,([Phi(i),X(1,i),X(2,i+1)]-[Phi(i),X(:,i)])/dt];
    %grad_p(:,i+1)=[(X(1,i+1)-X(1,i))/dt;(X(2,i+1)-X(2,i))/dt];
    %grad_Phi(:,i+1)=[(X(1,i+1)-X(1,i))/dt;(X(2,i+1)-X(2,i))/dt];
    grad_p(:,i+1)=[(p(i+1)-p(i))/Dxy(1,i);(p(i+1)-p(i))/Dxy(2,i)];
    grad_Phi(:,i+1)=[(Phi(i+1)-Phi(i))/Dxy(1,i);(Phi(i+1)-Phi(i))/Dxy(2,i)];
    grad_m(:,i+1)=[(m(1,i+1)-m(1,i))/Dxy(1,i);(m(2,i+1)-m(2,i))/Dxy(2,i)];
    laplace_p(i+1)=(grad_p(1,i+1)-grad_p(1,i))/Dxy(1,i)+(grad_p(2,i+1)-grad_p(2,i))/Dxy(2,i);
    laplace_Phi(i+1)=(grad_Phi(1,i+1)-grad_Phi(1,i))/Dxy(1,i)+(grad_Phi(2,i+1)-grad_Phi(2,i))/Dxy(2,i);
    %inv_neg_laplace_Phi(i+1)= ifft(-Phi(i+1))/(2*pi)^2*(X(1,i+1)^2+X(2,i+1)^2);%有问题
    inv_neg_laplace_Phi(i+1)= (-laplace_Phi(i+1))^-1;
    Fp = P_0*p(i+1)^3-3*P_0*(norm(m(:,i+1)-p(i+1)*A,2))^2*p(i+1)/U_tip^2-6*P_0*A'*(m(:,i+1)-p(i+1)*A)*p(i+1)^2/U_tip-d_0*p_a*s*B*(norm(m(:,i+1)-p(i+1)*A,2))^3-par_Phi_t*p(i+1)^3-eta_A*laplace_Phi(i+1)*p(i+1)^3+1/tau*(p(i+1)-p(i))*p(i+1)^3-1.5*d_0*p_a*s*B*norm(m(:,i+1)-p(i+1)*A,2)*A'*(m(:,i+1)-p(i+1)*A)*p(i+1); %（13）式的表达式
    DFp = 3*P_0*p(i+1)^2-3*P_0/U_tip^2*(-2*A'*(m(:,i+1)-p(i+1)*A)*p(i+1)+(norm(m(:,i+1)-p(i+1)*A,2))^2)-6*P_0*(-A'*A*p(i+1)^2+2*A'*(m(:,i+1)-p(i+1)*A)*p(i+1))/U_tip-d_0*p_a*s*B*((-2*A'*(m(:,i+1)-p(i+1)*A))*norm(m(:,i+1)-p(i+1)*A,2)+(norm(m(:,i+1)-p(i+1)*A,2))^2*(-A'*(m(:,i+1)-p(i+1)*A)/norm(m(:,i+1)-p(i+1)*A,2)))-3*par_Phi_t*p(i+1)^2-3*eta_A*laplace_Phi(i+1)*p(i+1)^2+(4*p(i+1)^3-3*p(i)*p(i+1)^2)/tau-1.5*d_0*p_a*s*B*(norm(m(:,i+1)-p(i+1)*A,2)*A'*(m(:,i+1)-p(i+1)*A)+p(i+1)*(-A'*(m(:,i+1)-p(i+1)*A)*A'*(m(:,i+1)-p(i+1)*A)/norm(m(:,i+1)-p(i+1)*A,2)-norm(m(:,i+1)-p(i+1)*A,2)*A'*A)); %（13）式的导函数表达式
    Fm = 6*P_0*(m(:,i+1)-p(i+1)*A)/(U_tip^2*p(i+1))+1.5*d_0*p_a*s*B*norm(m(:,i+1)-p(i+1)*A,2)*(m(:,i+1)-p(i+1)*A)/p(i+1)^2-grad_Phi(:,i+1)+1/tau*(m(:,i+1)-m(:,i));%14式子的表达式
    Gm= 6*P_0*p(i+1)/U_tip^2*eye(2)+1.5*d_0*p_a*s*B*((m(:,i+1)-p(i+1)*A)*(m(:,i+1)-p(i+1)*A)'/norm(m(:,i+1)-p(i+1)*A,2)+norm(m(:,i+1)-p(i+1)*A,2)*eye(2))+p(i+1)^2/tau*eye(2);%14式子的导数表达式
    %Phi(i+2) = Phi(i+1)+sigma*(par_p_t+divergence(m(:,i+1),X(:,i+1))-eta_A*laplace_p(i+1))*inv_neg_laplace_Phi(i+1); %%%/'Laplace算子';
    Phi(i+2) = Phi(i+1)+sigma*(par_p_t+grad_m(1,i+1)+grad_m(2,i+1)-eta_A*laplace_p(i+1))*inv_neg_laplace_Phi(i+1);
    
    p_star = p(i+1)- Fp/DFp;%牛顿迭代格式
    m_star= m(:,i+1)-[Gm(4),-Gm(3);-Gm(2),Gm(1)]/det(Gm)*Fm;
    Phi_star = Phi(i+2);
    

    if abs(p_star-p(i+1))<1e-4==1 && abs(m_star(1)-m(1,i+1))<1e-4==1 && abs(m_star(2)-m(2,i+1))<1e-4==1 && abs(Phi_star-Phi(i+1))<1e-4==1 %收敛判断，这里的阈值需要调试
        break;
    else
        p(i+2)=p_star;
        m(:,i+2)=m_star;
        Phi(i+2)=Phi_star;
    end
    i =i+1;

    fprintf('\n%s%.10f\t%s%d\n','p_star=',p_star,'i=',i);%输出p的结果
    fprintf('\n%s%.10f\t%s%d\n','m_star=',m_star,'i=',i);%输出m的结果
    fprintf('\n%s%.10f\t%s%d\n','Phi_star=',Phi_star,'i=',i);%输出Phi的结果
end
v_opt = m_star./p_star-A;
fprintf('\n%s%.10f\n','速度v_opt=',v_opt);
%fprintf('\n%s\n','Vxy=',Vxy);
disp(Vxy); 
%%作图 
 t=0.0001:dt:(I-1)*0.0001;
 plot(t,Vxy,'r');
 xlabel('时间');
 ylabel('UAV平均速度');
 legend('eta_a=1e-2');
 hold on