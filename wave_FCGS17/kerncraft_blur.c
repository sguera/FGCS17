float a[M][N];
float b[M][N];
float s0, s1, s2, s4, s5, s8, f;

for (int j = 2; j < M - 2; j++)
  for (int i = 2; i < N - 2; i++)
	b[j][i] = f * (
					s0 * a[j][i] +
					s1 * (a[j][i-1] + a[j][i+1] + a[j-1][i] + a[j+1][i]) +
					s2 * (a[j-1][i-1] + a[j-1][i+1] + a[j+1][i-1] + a[j+1][i+1]) +
					s4 * (a[j][i-2] + a[j][i+2] + a[j-2][i] + a[j+2][i]) +
					s5 * (a[j-1][i-2] + a[j-2][i-1] + a[j-2][i+1] + a[j-1][i+2] +
							a[j+1][i-2] + a[j+2][i-1] + a[j+2][i+1] + a[j+1][i+2]) +
					s8 * (a[j-2][i-2] + a[j-2][i+2] + a[j+2][i-2] + a[j+2][i+2])
				);

kerncraft -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml blur.c -p LC -p ECMData -p Roofline -D M 200 -D N 200 -D P 200
                                   kerncraft                                    
blur.c                                  -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml
-D P 200 -D M 200 -D N 200
-------------------------------------- LC --------------------------------------
1D layer condition:
L1: unconditionally fulfilled
L2: unconditionally fulfilled
L3: unconditionally fulfilled
2D layer condition:
L1: N <= 1366
L2: N <= 32770/3
L3: N <= 1572866/3

                                   kerncraft                                    
blur.c                                  -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml
-D P 200 -D M 200 -D N 200
----------------------------------- ECMData ------------------------------------
L1-L2 = 3.00 cy/CL
L2-L3 = 6.00 cy/CL
L3-MEM = 0.00 cy/CL

                                   kerncraft                                    
blur.c                                  -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml
-D P 200 -D M 200 -D N 200
----------------------------------- Roofline -----------------------------------
CPU bound with 1 cores(s)
21.36 GFLOP/s due to CPU max. FLOP/s