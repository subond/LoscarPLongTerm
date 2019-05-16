% Calculating statistics of Mrsay et al 2014 relationship between temperature and martin coefficient
T = [3.207396553
3.221575832
6.307072028
7.890334762
11.57119299
15.81139008
16.60615869
17.63488904]';
B = [0.491057125
0.570120601
0.880149191
0.685956661
0.700139084
1.35869387
1.244543396
1.591132307
]';
p = polyfit(T,B,1)
msre = sqrt(mean((B-(T.*0.062+0.303)).^2))
figure
scatter(T,B)
hold on
plot(T,T*p(1)+p(2))
plot(T, T*p(1)+p(2)+msre,'r')
plot(T, T*p(1)+p(2)-msre,'g')

