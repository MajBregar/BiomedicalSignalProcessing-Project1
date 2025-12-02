

% Basic Octave test script

% Generate data
t = 0:0.01:2*pi;
y = sin(t);

% Display a value
disp("Max value of sin(t):");
disp(max(y));

% Plot
figure;
plot(t, y);
title("Sine Wave");
xlabel("t");
ylabel("sin(t)");

