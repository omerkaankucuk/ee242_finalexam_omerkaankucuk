import numpy as np
import matplotlib.pyplot as plt
import mpcontrolPC

Ts = 0.0005
fs = 1/Ts
num = [2.866e+07]
den = [1, 8463, 4.592e+06]
gs = mpcontrolPC.tf(num,den)
gz = mpcontrolPC.c2d(gs,Ts)
print(gz)

[Kp, Ki, Kd] = mpcontrolPC.pid_tune(gs, "CHR", "PI")
print(Kp)
print(Ki)

P_s = mpcontrolPC.tf([Kp],[1])
P_z = mpcontrolPC.c2d(P_s,Ts)

I_s = mpcontrolPC.tf([Ki],[1,0])
I_z = mpcontrolPC.c2d(I_s,Ts)

C_z = mpcontrolPC.parallel(P_z, I_z)

print(C_z)

gol_z = mpcontrolPC.series(C_z,gz)
h = mpcontrolPC.tf([1],[1],Ts)
gcl_z = mpcontrolPC.feedback(gol_z,h)

#Estimated Output

N=4000
est_out = mpcontrolPC.step(gcl_z, 50, N)
t = np.arange(0,N/2000,1/2000)
plt.figure()
plt.plot(t,est_out)
plt.title("Estimated (Theoretical;) Output")
plt.xlabel("n")
plt.ylabel("Speed")
plt.show()

real_out = np.loadtxt("pi_stepout.dat")
plt.figure()
plt.plot(t,real_out)
plt.title("Real Output")
plt.xlabel("n")
plt.ylabel("Speed")
plt.show()

plt.figure()
plt.plot(t,est_out, label = "Estimated")
plt.plot(t,real_out, label = "Real")
plt.title("Real vs Estimated Output")
plt.xlabel("n")
plt.ylabel("Speed")
plt.legend()
plt.show()
