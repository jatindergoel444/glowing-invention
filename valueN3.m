clear all
clc
%-------------------------------------------------------------------------
Tt_jet=298;
Pt_jet=829000;
T_jet=248;
p_jet=438000;
u_jet=1200;
MW_jet=2;
R=8314;
rho_jet=p_jet*MW_jet/(R*T_jet);
Tt_air=1370;
Pt_air=333000;
T_air=1280;
p_air=261000;
u_air=458;
Y_O2=0.251; Y_N2=0.611; Y_H2O=0.138;
MW_air=1/((Y_O2/32)+(Y_N2/28)+(Y_H2O/18));
rho_air=p_air*MW_air/(R*T_air);
ru=1.96;
d_jet=2.49*10^(-3);
c1=1.6; c2=1/3; c3=1.3; c4=0.76; c5=0.009;
injection_length=0.358;
%-------------------------------------------------------------------------
x=0.37:0.001:0.5;
y=0:0.0001:0.0254;
z=-0.0190:0.001:0.0190;
clear x_cl
%x=0.359;

syms x_cl
tic
b=zeros;conc=zeros;f_cl=zeros;y_cl=zeros;n2=zeros;f=zeros;
for i=1:length(x)
    b(i)=d_jet*0.54*ru^(2/3)*((x(i)-injection_length)/d_jet)^c2;
    %b_neg(i)=d_jet*c5*ru^(2/3)*((x(i)-injection_length)/d_jet)^c2;
    conc(i)=c3*((rho_jet/rho_air)*(u_jet/u_air)^(-1)*((x(i)-injection_length)/d_jet)^(-2))^(1/3);
    rw=MW_jet/MW_air;
    f_cl(i)=conc(i)*rw/(1+(rw-1)*conc(i));
    y_cl(i)=d_jet*c4*((x(i)-injection_length)/d_jet)^c2*ru^(2/3);
    for j=1:length(y)
        parfor k=1:length(z)
            x1=x(i);y1=y(j);z1=z(k);
            dn=2*x_cl - 2*x1 - 2*c2*c4*ru^(2/3)*(y1 - c4*d_jet*ru^(2/3)*(-(injection_length - x_cl)/d_jet)^c2)*(-(injection_length - x_cl)/d_jet)^(c2 - 1);
            m(i,j,k)=double(solve(dn==0,x_cl));
            m_cl(i,j,k)=d_jet*c1*((m(i,j,k)-injection_length)/d_jet)^c2*ru^(2/3);
            n(i,j,k)=(x1-m(i,j,k))^2+(y1-m_cl(i,j,k))^2+z1^2;
            % syms x_cl
            % %n=(x-x_cl)^2+(y-c4*d_jet*ru^(2/3)*((x_cl-injection_length)/d_jet)^c2)^2+z^2;
            % dn=2*x_cl - 2*x1 - 2*c2*c4*ru^(2/3)*(y1 - c4*d_jet*ru^(2/3)*(-(injection_length - x_cl)/d_jet)^c2)*(-(injection_length - x_cl)/d_jet)^(c2 - 1);
            % %dn2=2*c2^2*c4^2*ru^(4/3)*(-(injection_length - x_cl)/d_jet)^(2*c2 - 2) - (2*c2*c4*ru^(2/3)*(y - c4*d_jet*ru^(2/3)*(-(injection_length - x_cl)/d_jet)^c2)*(c2 - 1)*(-(injection_length - x_cl)/d_jet)^(c2 - 2))/d_jet + 2;
            % m(i,j,k)=double(solve(dn==0,x_cl));
            % m_cl(i,j,k)=d_jet*c1*((m(i,j,k)-injection_length)/d_jet)^c2*ru^(2/3);
            % n(i,j,k)=sqrt((x1-m(i,j,k))^2+(y1-m_cl(i,j,k))^2+z1^2);
            n2(i,j,k)=(y1-y_cl(i))^2+z1^2;
            %error(i,j,k)=(abs(n2(i,j,k)-n(i,j,k))/n(i,j,k))*100;
            % if y1> y_cl(i,j,k)
            f(i,j,k)=f_cl(i)*exp(-n2(i,j,k)/(2*b(i)^2));
            f2(i,j,k)=f_cl(i)*exp(-n(i,j,k)/(2*b(i)^2));
           
            % else  
            %     f(i,j,k)=f_cl(i)*exp(-n(i,j,k)^2/(2*b_neg(i)^2));
            % end
        end
    end
end
toc
f_avg1=zeros;
for l=1:length(x)
    parfor m=1:length(y)
    f_plane=f(l,m,:);
    f_avg1(l,m,1)=mean(f_plane);
    end
end 
f_avg=mean(f_avg1,2);
parfor i=1:length(x)
df(i,:)=gradient(f_avg1(i,:),y);
end
toc
parfor s=1:length(x)
    mixture_var(s,:)=c5/c4*b(s).*abs(df(s,:));
end
toc