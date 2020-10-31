%
% Command file flysim.m for 3F1 Flight Control Experiment.
% Copyright: Cambridge University Engineering Department, October 1994.
% Author: M.C. Smith.
%



%*-------------------------2
% N = 10;
% M = 10;
% num=N;den=[1 M 0];	%  Numerator and denominator of plant 
% 			%  Laplace transfer function
% 
% runtime=5;   	% target simulation interval in seconds
% 
% wght=[5,0,0,0];	% entries are: impulse, step and sinusoid disturbance
		% weightings and sinusoidal frequency (Hz). Impulse and step
		% occur randomly between 0.2 and 0.6 secs. Sinusoid 
	 	% begins at t=0.
% Dtime = 0.3;
% Kgain = 2;
% num = Kgain*num;
% 
% P = tf(num,den,'InputDelay',0.3);
% h = bodeplot(P);
% nyd = nyquistplot(P);

%-------------------------------------2.1
% Dtime = 0.323;
% Kgain = 2;
% T = 4*Dtime/pi;
% num= sqrt(8)/Kgain;
% den=[T^3, 3*T^2, 3*T, 1];
% wght=[5*Kgain*Dtime,0,0,0];
% runtime=10;   	% target simulation interval in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2.2
% Dtime = 0.323;
% Kgain = 1;
% T = 4*Dtime/pi;
% num= num4;
% den=den4;
% wght=[0,0,1,0.66];
% runtime=10;   	% target simulation interval in seconds



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2.3
% Delay_add  = 0.323;
% K = 1;
% T = 0.45;
% num= 2*K;
% den=[T,-1];
% wght=[0.1,0,0,0];
% runtime=5;   	% target simulation interval in seconds
% sys = tf(num,den,'InputDelay',Delay_add);
% h = nyquistplot(sys);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2.4

% L = 0.5;
% g= 9.81;
% T = sqrt(L/g);
% 
% num= [T 1];
% den=[-T^2 0 1];
% sys = tf(num,den);
% h = nyquistplot(sys)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3
% num=[6.3,4.3,0.28];
% den=[1,11.2,19.6,16.2,0.91,0.27];
% wght=[10,0,0,0];
% runtime=60;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3
% 
% num=[6.3,4.3,0.28];
% den=[1,11.2,19.6,16.2,0.91,0.27];
% wght=[2,0,0,0];
% runtime=15;  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PIDDDDDDDD

% num=[6.3,4.3,0.28];
% den=[1,11.2,19.6,16.2,0.91,0.27];
% wght=[20,2,10,0.88];
% runtime=10;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samper=30;	% target sampling period in milliseconds

srate=(samper+1.3)/1000;	% anticipated average sampling period in secs
                            % was samper+0.6

grphc1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pid
% Kc = 17.4;
% Tc = 1.908;
% integ = 0;
% deriv = 0;
% yprev=0;
% Kp = 0.6*Kc;
% Ti=0.5*Tc;
% Td = 0.125*Tc;
for i=1:count
	set(hh,'Xdata',hx,'Ydata',hy+y*hz);
	pp=p(1,2);   %comment out for 3 
    %pp=-17.4*y;  %section3 
    
    %%%%%%%%%%%%pid*****
%     integ = -integ + y;
%     %integ=sign(integ)*min(abs(integ),2);
%     
%     deriv = (y-yprev)/srate;
%     pp= - Kp*(y+integ/Ti + Td*deriv);
%     yprev=y;
    
	pp=sign(pp)*min(max(0,abs(pp)-0.05),10);
    
    %3%%%%%%%%%%
    %pp=sign(pp)*min(max(0,abs(pp)-0),10);
    %3%%%%%%%%%%
    
	set(jh,'Xdata',jx,'Ydata',jy+pp*jz);
	drawnow;

	ylist(i)=y;
	ulist(i)=pp;

	x=adis*x + bdis*(pp+disturb(i));
	y=cdis*x + ddis*(pp+disturb(i));
	
	while (time2-time1<samper)
		time2=clock;time2=1000*(60*time2(5)+time2(6));
	end
	thetimes(i)=time2;
	time1=time2;
	
	if (y<-10 | y>10 )
		flg=1;crashind=i+1;
		thetimes(i+1)=thetimes(i)+samper;
		ylist(i+1)=y;
		ulist(i+1)=sign(p(1,2))*min(abs(p(1,2)),10);
		break;
	end

end

grphc2