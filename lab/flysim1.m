%
% Command file flysim.m for 3F1 Flight Control Experiment.
% Copyright: Cambridge University Engineering Department, October 1994.
% Author: M.C. Smith.
%



%Ask 
%-------------------------------------------------------------------------
%-------------------ASK FOR SECTION 2---------------------------------
%1------- 
%why in section 2, this will generate a different bode plot than bodedisp
% N = 10;
% M = 10;
% num=N;den=[1 M 0];	%  Numerator and denominator of plant 
% 
% % Dtime = 0.3;
% % Kgain = 2;
% % num = Kgain*num;
% % 
% % P = tf(num,den,'InputDelay',0.3);
% % h = bodeplot(P);
% % nyd = nyquistplot(P);
% P =
%  
%                     10
%   exp(-0.3*s) * ----------
%                 s^2 + 10 s
% 
%2------------
%Amount of extra time delay which can be tolerated 
% = Phase_margin(rad)/ws(freq when magnitude is zero) 
% = 37*(pi/180) /2 = 0.323 seconds  
%SEE https://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlFrequency
%  The phase margin also measures the system's tolerance to time delay. 
%  If there is a time delay greater than $PM/\omega_{gc}$ in the loop 
%  (where $\omega_{gc}$ is the frequency in rad/sec where the magnitude 
%  is 0 dB and PM is the phase margin converted to radians), 
%  the closed-loop system will become unstable. 
%  The time delay, $\tau_d$, can be thought of as an extra block 
%  in the forward path of the block diagram that adds phase lag to the system, 
%  but has no effect on the gain. 
%  That is, a time delay can be represented as a block with magnitude 
%  of 1 and phase $-\omega \tau_d$.




%-------------------ASK FOR SECTION 2.1---------------------------------
% How to read period of the zero input?
% 
% the period of oscillation with input is not constant? will the controller
% (my hand) affect its period?
% 
% Negative phase margin so unstable??
% (see: https://www.cds.caltech.edu/~murray/courses/cds101/fa02/faq/02-11-18_negativepm.html)

% how to find theoretical osc period? the freq at phase = -180??
% 
% the horizontal axis. Is it omega(=w*pi*f) ??

%guideline???? Avoid neg phase margin? Add a phase compensator???

%-------------------ASK FOR SECTION 2.2---------------------------------
% why can't I reduce the error. 
% Is my brain not powerful enough?(try this with pid)
% 
% Is the maximum propotional gain equals to gain margin?ASKKKKKKKKKKKKK
% MUST ASK 
% The gain margin is defined as the change in open-loop gain 
% required to make the closed-loop system unstable
% 
% 0.66hz -> 2pi*0.66, so find 4.146 on horizontal axis??

% Open loop tf is just G(s)?
% 
% Why move joystick more frequently?


%-------------------ASK FOR SECTION 2.3---------------------------------
% Nyquist diagram will have many curl at right half plane when delay added 
% but why?

% %-------------------ASK FOR SECTION 2.4---------------------------------
% Is my G(s) right. Why it is not possible to enclose -1
% L = 0.01;
% g= 9.81;
% T = sqrt(L/g);
% 
% num= [T 1];
% den=[-T^2 0 1];
% sys = tf(num,den);
% h = nyquistplot(sys)
%-----------------------------------------------------------------
% -------------------ASK FOR SECTION 3.2---------------------------------
% why can't I see the wind-up?
% Discuss the smallest number. 

% N = 10;
% M = 10;
% num=N;den=[1 M 0];	%  Numerator and denominator of plant 
% 			%  Laplace transfer function
% 
% runtime=5;   	% target simulation interval in seconds
% 
% wght=[5,0,0,0];	% entries are: impulse, step and sinusoid disturbance
% 		% weightings and sinusoidal frequency (Hz). Impulse and step
% 		% occur randomly between 0.2 and 0.6 secs. Sinusoid 
% 	 	% begins at t=0.



% Dtime = 0.03;
% Kgain = 3;
% T = 4*Dtime/pi;
% num= sqrt(8)/Kgain;
% den=[T^3, 3*T^2, 3*T, 1];
% wght=[5*Kgain*Dtime,0,0,0];
% runtime=10;   	% target simulation interval in seconds


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2.2
% Dtime = 0.03;
% Kgain = 1;
% T = 4*Dtime/pi;
% num= num4;
% den=den4;
% wght=[0,0,1,0.66];
% runtime=10;   	% target simulation interval in seconds



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2.3
% Delay_add  = 0.3;
% K = 0.8;
% T = 0.8;
% num= 2*K;
% den=[T,-1];
% wght=[0.1,0,0,0];
% runtime=5;   	% target simulation interval in seconds
% sys = tf(num,den,'InputDelay',Delay_add);
% h = nyquistplot(sys);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3
% num=[6.3,4.3,0.28];
% den=[1,11.2,19.6,16.2,0.91,0.27];
% wght=[10,0,0,0];
% runtime=60;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3

% num=[6.3,4.3,0.28];
% den=[1,11.2,19.6,16.2,0.91,0.27];
% wght=[2,0,0,0];
% runtime=15;  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PIDDDDDDDD
% 
num=[6.3,4.3,0.28];
den=[1,11.2,19.6,16.2,0.91,0.27];
wght=[0,2,0,0];
runtime=10;  



% num=num4;   %NOT THE RIGHT PID FOR THIS TF , SO PID NOT WORKING 
% den = den4;
% wght=[0,0,1,0.66];
% runtime=10;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samper=30;	% target sampling period in milliseconds

srate=(samper+1.3)/1000;	% anticipated average sampling period in secs
                            % was samper+0.6

grphc1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pid
Kc = 17.4;
Tc = 1.908;
integ = 0;
deriv = 0;
yprev=0;
Kp = 0.6*Kc;
Ti=0.5*Tc;
% Td = 0.125*Tc;
Td = 0.125*Tc*1.4;


for i=1:count
	set(hh,'Xdata',hx,'Ydata',hy+y*hz);

%     %%%%%%%%%%%%pid*****
    integ = integ + y*srate;
    integ=sign(integ)*min(abs(integ),0.175);
    %integ = 0;
    deriv = (y-yprev)/srate;
    pp= - Kp*(y+integ/Ti + Td*deriv);
    yprev=y;

    
    %%%%%%%%%%%%%%%manual control
%     pp=p(1,2);   
%     pp=sign(pp)*min(max(0,abs(pp)-0.05),10);

    
    
    
    
    
    
    
    
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