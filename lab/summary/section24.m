
L = 0.01;
g= 9.81;
T = sqrt(L/g);

num= [T 1];
den=[-T^2 0 1];
sys = tf(num,den);
h = nyquistplot(sys)