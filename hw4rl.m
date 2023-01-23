
K = 4.23;
sys = tf([0 0 0 1],[1 4 5 0]);

sys = tf([0 0 0 K],[1 4 5 K]);


stepplot(sys);