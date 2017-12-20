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