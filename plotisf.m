load('debug.txt');
load('a.txt');

load('dvdebug.txt');
load('dlambdadebug.txt');
load('dDdebug.txt');

tau=0.01:0.01:90;
figure(3);
plot((tau),debug);

figure(4);
plot(log10(abs(a)));

figure(5);
plot(log10(tau),dvdebug);

figure(6);
plot(log10(tau),dlambdadebug);

figure(7);
plot(log10(tau),dDdebug);