%%%%%%%%%%%%%%%%%%%%%%%%% Given information %%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Enter number of cycles for simulation'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
tt = str2double(answer(1)); 
dia = 20/1000; % diameter in m 
n = 4; % number of flutes
hang = deg2rad(30); % angle b/w flutes in radians
rang = deg2rad(10); % rake angle
sp = 2500/60; % rev/second
a = 1.5/1000; % axial depth of cut
fr = [0.025,0.050,0.100,0.150,0.200]/1000;
Fx = [-69.4665,-88.1012,-105.6237,-118.1107,-130.6807];
Fy = [64.3759,99.1383,155.7100,210.4644,263.4002];
Fz = [-10.2896,-20.5814,-33.9852,53.4833,73.1515];
c = fr; % given c = feed/rev 
k = 2*tan(hang)/dia;

%%%%%%%%%%%%%%%%%%%%%%%%% Finding cutting and edge constant %%%%%%%%%%%%%%%%%%%%%%%%%
cx = polyfit(fr,Fx,1);
cy = polyfit(fr,Fy,1);
cz = polyfit(fr,Fz,1);
x = linspace(0,0.225/1000);
y1 = cx(1)*(x) + cx(2);
y2 = cy(1)*(x) + cy(2);
y3 = cz(1)*(x) + cz(2);
figure('Name','Measured Data','NumberTitle','off');
hold on
plot(fr,Fx,'*r');
plot(x,y1,'r');
plot(fr,Fy,'*g');
plot(x,y2,'g');
plot(fr,Fz,'*b');
plot(x,y3,'b');
legend('Fx data','Fx','Fy data','Fy','Fz data','Fz')
xlabel('feed rate (mm/rev)') 
ylabel('Force(N)') 
hold off
Krc = -cx(1)*4/(n*a); Kre = -cx(2)*pi/(n*a);
Ktc = cy(1)*4/(n*a); Kte = cy(2)*pi/(n*a);
Kac = -cz(1)*pi/(n*a); Kae = cz(2)*2/(n*a);

%%%%%%%%%%%%%%%%%%%%%%%%% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%
theta1 = linspace(0,pi);
zat1 = linspace(0,a,3);
meanx = []; meany = []; meanz = [];
for i = 1:length(c)
    Fx_og = @(theta,zat) (0.5*c(i)*(-Ktc*sin(2*theta-2*k*zat)-Krc*(1-cos(2*theta-2*k*zat))) + (-Kte*cos(theta-k*zat)-Kre*sin(theta-k*zat)) );
    Fy_og = @(theta,zat) (0.5*c(i)*(Ktc*(1-cos(2*theta-2*k*zat))-Krc*sin(2*theta-2*k*zat)) + (Kte*sin(theta-k*zat)-Kre*cos(theta-k*zat)) );
    Fz_og = @(theta,zat) (Kac*c(i)*sin(theta-k*zat)+Kae);
    
    meanx = [meanx ; integral2(Fx_og,0,pi,0,a)*n/(2*pi)];
    meany = [meany ; integral2(Fy_og,0,pi,0,a)*n/(2*pi)];
    meanz = [meanz ; integral2(Fz_og,0,pi,0,a)*n/(2*pi)];

    figure('Name',strcat('Question 2 cutting_differential_data_',string(i)),'NumberTitle','off')
    subplot(1,3,1)
    hold on
    plot(theta1,Fx_og(theta1,zat1(1)),'r')
    plot(theta1,Fx_og(theta1,zat1(2)),'g')
    plot(theta1,Fx_og(theta1,zat1(3)),'b')
    legend('z = 0','z = a/2','z = a')
    xlabel('theta') 
    ylabel('Force(N)') 
    title('differential Fx')

    subplot(1,3,2)
    hold on
    plot(theta1,Fy_og(theta1,zat1(1)),'r')
    plot(theta1,Fy_og(theta1,zat1(2)),'g')
    plot(theta1,Fy_og(theta1,zat1(3)),'b')
    legend('z = 0','z = a/2','z = a')
    xlabel('theta(in radians)') 
    ylabel('Force(N)') 
    title('differential Fy')
    
    subplot(1,3,3)
    hold on
    plot(theta1,Fz_og(theta1,zat1(1)),'r')
    plot(theta1,Fz_og(theta1,zat1(2)),'g')
    plot(theta1,Fz_og(theta1,zat1(3)),'b')
    legend('z = 0','z = a/2','z = a')
    xlabel('theta') 
    ylabel('Force(N)') 
    title('differential Fz')
end
figure('Name','Results','NumberTitle','off');
hold on
plot(fr,Fx,'*r');
plot(x,y1,'r');
plot(fr,meanx,'or');
plot(fr,Fy,'*g');
plot(x,y2,'g');
plot(fr,meany,'og');
plot(fr,Fz,'*b');
plot(x,y3,'b');
plot(fr,meanz,'ob');
legend('Fx measured','Fx','Fx theoretical','Fy measured','Fy','Fy theoretical','Fz measured','Fz','Fz theoretical')
xlabel('feed rate (mm/rev)') 
ylabel('Force(N)') 
hold off

%%%%%%%%%%%%%%%%%%%%%%%%% Question 3 part(i) %%%%%%%%%%%%%%%%%%%%%%%%%
meanxf1 =[];meanyf1 = []; meanzf1 = [];
phi = [];
for i = 1:(tt*4)
    phi = [phi, linspace(0,pi/2)];
end
for i = 1:length(c)
    Fx_og = @(theta,zat) (0.5*c(i)*(-Ktc*sin(2*theta-2*k*zat)-Krc*(1-cos(2*theta-2*k*zat))) + (-Kte*cos(theta-k*zat)-Kre*sin(theta-k*zat)) );
    Fy_og = @(theta,zat) (0.5*c(i)*(Ktc*(1-cos(2*theta-2*k*zat))-Krc*sin(2*theta-2*k*zat)) + (Kte*sin(theta-k*zat)-Kre*cos(theta-k*zat)) );
    Fz_og = @(theta,zat) (Kac*c(i)*sin(theta-k*zat)+Kae);
    
    meanxf1 = [meanxf1 ; integral2(Fx_og,0,pi/2,0,a)*n/(2*pi)];
    meanyf1 = [meanyf1 ; integral2(Fy_og,0,pi/2,0,a)*n/(2*pi)];
    meanzf1 = [meanzf1 ; integral2(Fz_og,0,pi/2,0,a)*n/(2*pi)];

    Fx_int = []; Fy_int = []; Fz_int = [];
    for cnt = 1:length(phi) 
        Fx_og1 = @(it) Fx_og(phi(cnt),it);
        Fy_og1 = @(it) Fy_og(phi(cnt),it);
        Fz_og1 = @(it) Fz_og(phi(cnt),it);
        Fx_int(cnt) = integral(Fx_og1,0,a);
        Fy_int(cnt) = integral(Fy_og1,0,a);
        Fz_int(cnt) = integral(Fz_og1,0,a);
    end
    figure('Name','simulated cycle motion','NumberTitle','off');
    hold on
    plot((2*pi)/400*linspace(1,tt*400,tt*400),Fx_int,'r')
    plot((2*pi)/400*linspace(1,tt*400,tt*400),Fy_int,'g')
    plot((2*pi)/400*linspace(1,tt*400,tt*400),Fz_int,'b')
    legend('Fx','Fy','Fz')
    xlabel('theta(in radians)') 
    ylabel('Force(N)') 
    title(strcat('Final forces for c = ',string(c(i)*1000)))
    hold off
end
cxf = polyfit(fr,meanxf1,1);
cyf = polyfit(fr,meanyf1,1);
czf = polyfit(fr,meanzf1,1);
y1f = cxf(1)*(x) + cxf(2);
y2f = cyf(1)*(x) + cyf(2);
y3f = czf(1)*(x) + czf(2);
figure('Name','Results_3(i)','NumberTitle','off');
hold on
plot(fr,meanxf1,'*r');
plot(x,y1f,'r--');
plot(x,y1,'r');
plot(fr,meanx,'or');

plot(fr,meanyf1,'*g');
plot(x,y2f,'g--');
plot(x,y2,'g');
plot(fr,meany,'og');

plot(fr,meanzf1,'*b');
plot(x,y3f,'b--');
plot(x,y3,'b');
plot(fr,meanz,'ob');

legend('Fx simulated data','Fx simulated','Fx','Fx theoretical','Fy simulated data','Fy simulated','Fy','Fy theoretical','Fz simulated data','Fz simulated','Fz','Fz theoretical')
xlabel('feed rate (mm/rev)') 
ylabel('Force(N)') 
hold off

%%%%%%%%%%%%%%%%%%%%%%%%% Question 3 part(ii) %%%%%%%%%%%%%%%%%%%%%%%%%
meanxf2 =[];meanyf2 = []; meanzf2 = [];
phi = [];
for i = 1:(tt*4)
    phi = [phi, linspace(pi/2,0)];
end
for i = 1:length(c)
    Fx_og = @(theta,zat) (0.5*c(i)*(-Ktc*sin(2*theta-2*k*zat)-Krc*(1-cos(2*theta-2*k*zat))) + (-Kte*cos(theta-k*zat)-Kre*sin(theta-k*zat)) );
    Fy_og = @(theta,zat) (0.5*c(i)*(Ktc*(1-cos(2*theta-2*k*zat))-Krc*sin(2*theta-2*k*zat)) + (Kte*sin(theta-k*zat)-Kre*cos(theta-k*zat)) );
    Fz_og = @(theta,zat) (Kac*c(i)*sin(theta-k*zat)+Kae);
    
    meanxf2 = [meanxf2 ; integral2(Fx_og,pi/2,0,0,a)*n/(2*pi)];
    meanyf2 = [meanyf2 ; integral2(Fy_og,pi/2,0,0,a)*n/(2*pi)];
    meanzf2 = [meanzf2 ; integral2(Fz_og,pi/2,0,0,a)*n/(2*pi)];

    Fx_int = []; Fy_int = []; Fz_int = [];
    for cnt = 1:length(phi) 
        Fx_og1 = @(it) Fx_og(phi(cnt),it);
        Fy_og1 = @(it) Fy_og(phi(cnt),it);
        Fz_og1 = @(it) Fz_og(phi(cnt),it);
        Fx_int(cnt) = integral(Fx_og1,0,a);
        Fy_int(cnt) = integral(Fy_og1,0,a);
        Fz_int(cnt) = integral(Fz_og1,0,a);
    end
    figure('Name','simulated cycle motion (ii)','NumberTitle','off');
    hold on
    plot((2*pi)/400*linspace(1,tt*400,tt*400),Fx_int,'r')
    plot((2*pi)/400*linspace(1,tt*400,tt*400),Fy_int,'g')
    plot((2*pi)/400*linspace(1,tt*400,tt*400),Fz_int,'b')
    legend('Fx','Fy','Fz')
    xlabel('theta(in radians)') 
    ylabel('Force(N)') 
    title(strcat('Final forces for c = ',string(c(i)*1000)))
    hold off
end
cxf = polyfit(fr,meanxf2,1);
cyf = polyfit(fr,meanyf2,1);
czf = polyfit(fr,meanzf2,1);
y1f = cxf(1)*(x) + cxf(2);
y2f = cyf(1)*(x) + cyf(2);
y3f = czf(1)*(x) + czf(2);
figure('Name','Results_3(ii)','NumberTitle','off');
hold on
plot(fr,meanxf2,'*r');
plot(x,y1f,'r--');
plot(x,y1,'r');
plot(fr,meanx,'or');

plot(fr,meanyf2,'*g');
plot(x,y2f,'g--');
plot(x,y2,'g');
plot(fr,meany,'og');

plot(fr,meanzf2,'*b');
plot(x,y3f,'b--');
plot(x,y3,'b');
plot(fr,meanz,'ob');

legend('Fx simulated data','Fx simulated','Fx','Fx theoretical','Fy simulated data','Fy simulated','Fy','Fy theoretical','Fz simulated data','Fz simulated','Fz','Fz theoretical')
xlabel('feed rate (mm/rev)') 
ylabel('Force(N)') 
hold off
