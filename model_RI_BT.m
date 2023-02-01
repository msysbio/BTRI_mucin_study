function model_RI_BT

%%% Author: Didier Gonze

%%% Simulations in monoculture (BT alone + RI alone) and in coculture (BT+RI)
%%% (1) in absence of mucin (mucin=0) => figure(1)
%%% (2) in presence of mucin (mucin=1) => figure(2)


clear
clc

for mucin=0:1

S0=mucin;   % mucin (0 = absent or 1 = present)


%%% Initial conditions

X1=0.01;  % RI
X1m=0;    % RI (slow-growth)
X2=0.01;  % BT
X2m=0;    % BT (slow-growth)

S1=1;     % glucose
S2=0;     % acetate/lactate
S3=0.1;   % mucin sugars
S4=0;     % butyrate

ICs=[0  0   X2 X2m S1 S2 S3 S4;    % BT mono-culture
     X1 X1m 0  0 S1 S2 S3 S4;      % RI mono-culture
     X1 X1m X2 X2m S1 S2 S3 S4];   % BT+RI co-culture


%%% Time

tend=120;
tstep=0.1;


%%% Figure

figure(mucin+1)
set(figure(mucin+1),'position',[550 80 1100 900])
clf


%%% Loop on IC

for k=1:3

IC=ICs(k,:);


%%% Integration

[t v] = ode45(@dxdt,[0:tstep:tend],IC,[],S0);

X1=v(:,1);   % RI (growth)
X1m=v(:,2);  % RI (slow-growth)
X2=v(:,3);   % BT
X2m=v(:,4);  % BT (slow-growth)

S1=v(:,5);  % glucose
S2=v(:,6);  % acetate/lactate
S3=v(:,7);  % mucin sugars
S4=v(:,8);  % butyrate


%%% pH functions

P=fph(S2);
phi1=P(:,1);
phi2=P(:,2);


%%% Growth functions

mu11=fmu11(S1,S3);
mu13=fmu13(S3);
mum12 = fmum12(S2);

mu1=mu11+mu13;

mu21 = fmu21(S1);
mu23 = fmu23(S3,S1);

mu2=mu21+mu23;


%%% Figure

subplot(4,3,k)

if (k==2 || k==3); plot(t,X1,'b','linewidth',1,'color',[0.8 0 0.5]);
    hold on; plot(t,X1m,'b--','linewidth',1,'color',[0.8 0 0.5]);
    hold on; plot(t,X1+X1m,'b','linewidth',2,'color',[0.8 0 0.5]);
end
hold on
if (k==1 || k==3); plot(t,X2,'r','linewidth',1,'color',[0 0.5 0.3]); 
    hold on; plot(t,X2m,'r--','linewidth',1,'color',[0 0.5 0.3]);
    hold on; plot(t,X2+X2m,'r','linewidth',2,'color',[0 0.5 0.3]);
end
if (k==3); plot(t,X1+X1m+X2+X2m,'k','linewidth',1,'color',[0.6 0.6 0.6]); 
end

ylim([0 0.22]);
xlim([0 tend])
set(gca,'xtick',[0:20:tend],'fontsize',14);
xlabel('Time','fontsize',16)
ylabel('X1, X2','fontsize',16)
if (k==1); title('Monoculture BT','fontsize',18); 
elseif (k==2); title('Monoculture RI','fontsize',18);
else (k==3); title('Coculture BT+RI','fontsize',18); 
L=legend('RI (growing)','RI (slow-growth)','RI (total)','BT (growing)','BT (slow-growth)','BT (total)','RI + BT (total)');
set(L,'fontsize',11);
end
box on;

subplot(4,3,3+k)

plot([0 tend],[S0 S0],'k','linewidth',2)
hold on
plot(t,S1,'linewidth',2,'color','b')
hold on
plot(t,[S2 S3],'linewidth',2)
hold on
plot(t,S4,'r','linewidth',2,'color',[0.9 0.8 0.8])

ylim([0 1.3]);
xlim([0 tend])
set(gca,'xtick',[0:20:tend],'fontsize',14);
xlabel('Time','fontsize',16)
ylabel('S0, S1, S2, S3, S4','fontsize',16)
if k==3    
L=legend('S0 = mucin', 'S1 = glucose','S2 = acetate/lactate','S3 = mucin sugars','S4 = butyrate');
set(L,'fontsize',11);
end
box on;


subplot(4,3,6+k)

zz=zeros(length(t),1);
if k==1
plot(t,[zz zz zz mu21 mu23],'linewidth',2)
elseif k==2
plot(t,[mu11 mu13 mum12 zz zz],'linewidth',2)
else
plot(t,[mu11 mu13 mum12 mu21 mu23],'linewidth',2)
end

xlim([0 tend])
set(gca,'xtick',[0:20:tend],'fontsize',14);
xlabel('Time','fontsize',16)
ylabel('Growth rates','fontsize',16)
if k==3
L=legend('\mu11','\mu13','\mu12','\mu21','\mu23');
set(L,'fontsize',10);
end
box on;


subplot(4,3,9+k)

if (k==2 || k==3); plot(t,phi1,'b','linewidth',2); end
hold on
if (k==1 || k==3); plot(t,phi2,'r','linewidth',2); end
ylim([0 1.1]);
xlim([0 tend])
set(gca,'xtick',[0:20:tend],'fontsize',14);
xlabel('Time','fontsize',16)
ylabel('\phi1, \phi2','fontsize',16)
if k==3
L=legend('\phi1','\phi2');
set(L,'fontsize',12);
end
box on;


%%% Save data

R=[t v];

outfile = sprintf('results_S0_%d_k_%d.dat',S0,k);
save(outfile,'R','-ASCII');

fprintf('results saved in %s ',outfile)
if (k==1) fprintf('(BT monoculture, ')
elseif (k==2) fprintf('(RI monoculture, ')
else fprintf('(RI/BT coculture, ')
end 
if (S0==0) fprintf('no mucin)')
elseif (S0==1) fprintf('with mucin)')
end    
fprintf('\n',outfile)



end  % end loop on k

end  % end loop on mucin



% ======================================================================
% Growth functions for RI
% ======================================================================

function mu11 = fmu11(S1,S3)

v11=0.6;   K11=0.5;       % growth of RI on glucose
ki11=1;                   % inhibition by mucin sugars of RI growth on glucose

finhib1=ki11./(ki11+S3);

mu11=v11*finhib1.*S1./(K11+S1);


function mu13 = fmu13(S3)

v13=0.1;  K13=0.5;        % growth of RI on other mucin sugars

mu13=v13*S3./(K13+S3);


function mum12 = fmum12(S2)

v12=0.01; K12=0.1;        % slow growth of RI on acetate/lactate

mum12=v12.*S2./(K12+S2);  % growth rate of RI on acetate/lactate


% ======================================================================
% Growth functions for BT
% ======================================================================

function mu21 = fmu21(S1)

v21=0.7;   K21=0.5;       % growth of BT on glucose

mu21=v21*S1./(K21+S1);


function mu23 = fmu23(S3,S1)

v23=0.6;   K23=0.05;      % growth of BT on other mucin sugars
Ki13=0.1;                 % inhibition by glucose of BT growth on mucin sugars  
n=4;

finhib=Ki13.^n./(Ki13.^n+S1.^n);  

mu23=finhib.*v23.*S3./(K23+S3);


% ======================================================================
% Switch functions for RI
% ======================================================================

function km1 = fkm1(S1,S3)

vm1=0.04;            % rate of switch from growing RI to slow-growth mode
Ki1=0.1; Ki3=0.6;    % threshold <-> switch from growing RI to slow-growth mode
h=4;                 % Hill coefficient <-> switch of RI to slow-growth mode

km1=vm1*(Ki1^h/(Ki1^h+S1^h))*(Ki3^h/(Ki3^h+S3^h));  % rate of switch to slow-growth state for RI (depends on glucose and mucin sugars)


% ======================================================================
% Switch functions for BT
% ======================================================================

function km2 = fkm2(S0)

vm20=0.001;       % basal rate of switch of BT to slow-growth mode
vm2S=0.1;         % rate of switch of BT to slow-growth mode boosted by mucin

km2=vm20+vm2S*S0; % rate of switch to slow-growth state for BT (boosted by mucin)


% ======================================================================
% pH functions
% ======================================================================

function P = fph(S2)

%%% S2 = acid (lactate/acetate)

ki1=0.8;    % "sensitivity" of RI to pH
ki2=0.6;    % "sensitivity" of BT to pH
m=2;

phi1=ki1.^m./(ki1.^m+S2.^m);  % for RI
phi2=ki2.^m./(ki2.^m+S2.^m);  % for BT

P=[phi1 phi2];


% ======================================================================
% Equations
% ======================================================================

function dv = dxdt(t,v,S0)

%%% Variables

X1=v(1);    % RI 
X1m=v(2);   % RI (slow-growth)
X2=v(3);    % BT
X2m=v(4);   % BT (slow-growth)

S1=v(5);    % glucose
S2=v(6);    % acetate/lactate
S3=v(7);    % mucin sugars
S4=v(8);    % butyrate


%%% Parameters

d1=0.2;        % death rate of RI
d2=0.2;        % death rate of BT

f1=0.25;      % fraction of attached RI
f2=0.25;      % fraction of attached BT

dm1=0.004;    % death rate of RI in slow-growth state (very low) 
dm2=0.005;    % death rate of BT in slow-growth state (very low) 

g11=2; g13=2; g12m=2;  % consumption rates (<-> yields for RI)
g21=2; g23=0.5;        % consumption rates (<-> yields for BT)

a12=1; a13=2; a14=2; a14m=2;    % production rates (<-> RI)
a22=1; a23=2;                   % production rates (<-> BT)


%%% attached fractions

Xb1=f1*S0*X1;
Xb2=f2*S0*X2;


%%% pH functions

P=fph(S2);
phi1=P(:,1);  % for RI
phi2=P(:,2);  % for BT


%%% growth functions

mu11=fmu11(S1,S3);
mu13=fmu13(S3);

mu21=fmu21(S1);
mu23=fmu23(S3,S1);

mum12=fmum12(S2);

mu1=phi1*(mu11+mu13);
mu2=phi2*(mu21+mu23);

mum1=phi1*mum12;
mum2=0; 


%%% switches to slow-growth modes

km1 = fkm1(S1,S3);
km2 = fkm2(S0);


%%% Equations

dv = [
    (mu1-d1-km1)*X1;                    % dX1/dt (RI) 
    km1*X1+(mum12-dm1)*X1m;             % dX1m/dt (RI slow-growth)
    (mu2-d2-km2)*X2;                    % dX2/dt (BT)
    km2*X2+(mum2-dm2)*X2m;              % dX2m/dt (BT slow-growth)
    -g11*mu11*X1-g21*mu21*X2;                       % dS1/dt (glucose)
    a12*mu11*X1+a22*mu21*X2-g12m*mum12*X1m;         % dS2/dt (acetate/lactate)
    a13*S0*Xb1+a23*S0*Xb2-g13*mu13*X1-g23*mu23*X2;  % dS3/dt (mucin sugars)
    a14*mu1*X1+a14m*mum1*X1m;                       % dS4/dt (butyrate)
] ; 






