clear all
close all
clc

[TTS2] = funExercise2(120*ones(1,60));
[X31,Fval31,flag31] = fmincon(@funExercise2,120*ones(1,60),[],[],[],[],60*ones(1,60),120*ones(1,60));
[X32,Fval32,flag32] = fmincon(@funExercise2,60*ones(1,60),[],[],[],[],60*ones(1,60),120*ones(1,60));
[X4,Fval4,flag4] = ga(@funExercise4, 120, [],[],[],[],[60*ones(1,60) zeros(1,60)],[120*ones(1,60), ones(1,60)]);
[X6,Fval6,flag6] = ga(@funExercise6, 120, [],[],[],[],[60*ones(1,60) zeros(1,60)],[120*ones(1,60), ones(1,60)]);

[k2,wr2,r2,Vsl2] = dataVsl(120*ones(1,60));
[k31,wr31,r31,Vsl31] = dataVsl(X31);
[k32,wr32,r32,Vsl32] = dataVsl(X32);
[k4,wr4,r4,Vsl4] = dataVslAndR(X4);

T = 10;
figure
plot((0:k2)*T,wr2,'-',(0:k31)*T,wr31,'-',(0:k32)*T,wr32,'-',(0:k4)*T,wr4,'-')
title('queue on ramp vs time')
xlabel('Time [s]')
ylabel('Cars in queue [cars]')
legend('uncontrolled','Vsl controlled, initial 120', 'Vsl controlled, initial 60', 'Vsl and r controlled')
figure
plot((1:k2)*T,r2,'-',(1:k31)*T,r31,'-',(1:k32)*T,r32,'-',(1:k4)*T,r4,'-')
title('ramp metering rate vs time')
xlabel('Time [s]')
ylabel('Ramp metering rate')
legend('uncontrolled','Vsl controlled, initial 120', 'Vsl controlled, initial 60', 'Vsl and r controlled')
figure
plot((1:k2)*T,Vsl2,'-',(1:k31)*T,Vsl31,'-',(1:k32)*T,Vsl32,'-',(1:k4)*T,Vsl4,'-')
title('Vsl vs time')
xlabel('Time [s]')
ylabel('Vsl [km/h]')
legend('uncontrolled','Vsl controlled, initial 120', 'Vsl controlled, initial 60', 'Vsl and r controlled')

function[TTS] = funExercise2(Vsl)
r=ones(1,60);

E1 = 3;
E2 = 0.5;
E3 = 4.5;

tau = 10;
mu = 80;
Cr = 2000;
rhom = 120;
alfa = 0.1;
K = 10;
a = 2;
vf = 110;
rhoc = 28;
T = 10/3600;
lambda = 3;
L = 1;


xr(:,1) = [20; 20; 20; 20];
xv(:,1) = [90; 90; 90; 90];
wr(1) = 0;

for k=1:60
    
    Dr(k) = 1500;
    qr(k) = min([r(k)*Cr Dr(k)+wr(k)/T Cr*(rhom - xr(4,k))/(rhom - rhoc)]);
    if k<12
        q0(k) = (7000 + 100*E1);
    else
        q0(k) = (2000 + 100*E2);
    end
    
    for i=1:4
        q(i,k) = lambda * xr(i,k) * xv(i,k);
        if i == 1
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q0(k) - q(i,k));
        elseif i == 4
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k) + qr(k));
        else
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k));
        end
    end
    
    for i=1:4
        if i == 1
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        elseif  i == 4
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i,k)-xr(i,k))/(xr(i,k)+K));
        else
            V(i,k) = min([(1+alfa)*Vsl(k), vf*exp((-1/a)*(xr(i,k)/rhoc)^a)]);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        end
    end
    
    wr(k+1) = wr(k) + T*(Dr(k)-qr(k));
    y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k));
end
TTS = sum(y);
end

function[k,wr,r,Vsl] = dataVsl(Vsl)
r=ones(1,60);

E1 = 3;
E2 = 0.5;
E3 = 4.5;

tau = 10;
mu = 80;
Cr = 2000;
rhom = 120;
alfa = 0.1;
K = 10;
a = 2;
vf = 110;
rhoc = 28;
T = 10/3600;
lambda = 3;
L = 1;


xr(:,1) = [20; 20; 20; 20];
xv(:,1) = [90; 90; 90; 90];
wr(1) = 0;

for k=1:60
    
    Dr(k) = 1500;
    qr(k) = min([r(k)*Cr Dr(k)+wr(k)/T Cr*(rhom - xr(4,k))/(rhom - rhoc)]);
    if k<12
        q0(k) = (7000 + 100*E1);
    else
        q0(k) = (2000 + 100*E2);
    end
    
    for i=1:4
        q(i,k) = lambda * xr(i,k) * xv(i,k);
        if i == 1
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q0(k) - q(i,k));
        elseif i == 4
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k) + qr(k));
        else
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k));
        end
    end
    
    for i=1:4
        if i == 1
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        elseif  i == 4
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i,k)-xr(i,k))/(xr(i,k)+K));
        else
            V(i,k) = min([(1+alfa)*Vsl(k), vf*exp((-1/a)*(xr(i,k)/rhoc)^a)]);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        end
    end
    
    wr(k+1) = wr(k) + T*(Dr(k)-qr(k));
    if wr(k) > 20-E3
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k))+1e10;
    else
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k));
    end
end
TTS = sum(y);
end

function[TTS] = funExercise4(X)
Vsl = X(1:60);
r = X(61:120);

E1 = 3;
E2 = 0.5;
E3 = 4.5;

tau = 10;
mu = 80;
Cr = 2000;
rhom = 120;
alfa = 0.1;
K = 10;
a = 2;
vf = 110;
rhoc = 28;
T = 10/3600;
lambda = 3;
L = 1;


xr(:,1) = [20; 20; 20; 20];
xv(:,1) = [90; 90; 90; 90];
wr(1) = 0;

for k=1:60
    
    Dr(k) = 1500;
    qr(k) = min([r(k)*Cr Dr(k)+wr(k)/T Cr*(rhom - xr(4,k))/(rhom - rhoc)]);
    if k<12
        q0(k) = (7000 + 100*E1);
    else
        q0(k) = (2000 + 100*E2);
    end
    
    for i=1:4
        q(i,k) = lambda * xr(i,k) * xv(i,k);
        if i == 1
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q0(k) - q(i,k));
        elseif i == 4
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k) + qr(k));
        else
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k));
        end
    end
    
    for i=1:4
        if i == 1
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        elseif  i == 4
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i,k)-xr(i,k))/(xr(i,k)+K));
        else
            V(i,k) = min([(1+alfa)*Vsl(k), vf*exp((-1/a)*(xr(i,k)/rhoc)^a)]);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        end
    end
    
    wr(k+1) = wr(k) + T*(Dr(k)-qr(k));
    if wr(k) > 20-E3
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k))+1e10;
    else
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k));
    end
end

TTS = sum(y);
end

function[k,wr,r,Vsl] = dataVslAndR(X)
Vsl = X(1:60);
r = X(61:120);

E1 = 3;
E2 = 0.5;
E3 = 4.5;

tau = 10;
mu = 80;
Cr = 2000;
rhom = 120;
alfa = 0.1;
K = 10;
a = 2;
vf = 110;
rhoc = 28;
T = 10/3600;
lambda = 3;
L = 1;


xr(:,1) = [20; 20; 20; 20];
xv(:,1) = [90; 90; 90; 90];
wr(1) = 0;

for k=1:60
    
    Dr(k) = 1500;
    qr(k) = min([r(k)*Cr Dr(k)+wr(k)/T Cr*(rhom - xr(4,k))/(rhom - rhoc)]);
    if k<12
        q0(k) = (7000 + 100*E1);
    else
        q0(k) = (2000 + 100*E2);
    end
    
    for i=1:4
        q(i,k) = lambda * xr(i,k) * xv(i,k);
        if i == 1
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q0(k) - q(i,k));
        elseif i == 4
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k) + qr(k));
        else
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k));
        end
    end
    
    for i=1:4
        if i == 1
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        elseif  i == 4
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i,k)-xr(i,k))/(xr(i,k)+K));
        else
            V(i,k) = min([(1+alfa)*Vsl(k), vf*exp((-1/a)*(xr(i,k)/rhoc)^a)]);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        end
    end
    
    wr(k+1) = wr(k) + T*(Dr(k)-qr(k));
    if wr(k) > 20-E3
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k))+1e10;
    else
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k));
    end
end

TTS = sum(y);
end

function[TTS] = funExercise6(X)
Vsl = X(1:60);

for j = 1:60
    if Vsl(j) < 70
        Vsld(j) = 60;
    elseif Vsl(j) < 90
        Vsld(j) = 80;
    elseif Vsl(j) < 110
        Vsld(j) = 100;
    else
        Vsld(j) = 120;
    end
end
        
r = X(61:120);

E1 = 3;
E2 = 0.5;
E3 = 4.5;

tau = 10;
mu = 80;
Cr = 2000;
rhom = 120;
alfa = 0.1;
K = 10;
a = 2;
vf = 110;
rhoc = 28;
T = 10/3600;
lambda = 3;
L = 1;


xr(:,1) = [20; 20; 20; 20];
xv(:,1) = [90; 90; 90; 90];
wr(1) = 0;

for k=1:60
    
    Dr(k) = 1500;
    qr(k) = min([r(k)*Cr Dr(k)+wr(k)/T Cr*(rhom - xr(4,k))/(rhom - rhoc)]);
    if k<12
        q0(k) = (7000 + 100*E1);
    else
        q0(k) = (2000 + 100*E2);
    end
    
    for i=1:4
        q(i,k) = lambda * xr(i,k) * xv(i,k);
        if i == 1
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q0(k) - q(i,k));
        elseif i == 4
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k) + qr(k));
        else
            xr(i,k+1) = xr(i,k) + (T/(lambda*L)) * (q(i-1,k) - q(i,k));
        end
    end
    
    for i=1:4
        if i == 1
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        elseif  i == 4
            V(i,k) = vf*exp((-1/a)*(xr(i,k)/rhoc)^a);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i,k)-xr(i,k))/(xr(i,k)+K));
        else
            V(i,k) = min([(1+alfa)*Vsld(k), vf*exp((-1/a)*(xr(i,k)/rhoc)^a)]);
            xv(i,k+1) = xv(i,k) + (T/tau) * (V(i,k)-xv(i,k)) + (T/L)*xv(i,k) * (xv(i-1,k)-xv(i,k)) - ((mu*T)/(tau*L))*((xr(i+1,k)-xr(i,k))/(xr(i,k)+K));
        end
    end
    
    wr(k+1) = wr(k) + T*(Dr(k)-qr(k));
    if wr(k) > 20-E3
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k))+1e10;
    else
        y(k) = T*wr(k) + T*L*lambda*sum(xr(:,k));
    end
end

TTS = sum(y);
end