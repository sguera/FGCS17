void wave_parm(float *  *  U_0_1_out, float *  U_0_m1, float *  U_0_0, float *  U_0_1, float fMin, float fDX, float fDT_DX_sq, long X_MAX, long Y_MAX, long Z_MAX, long cb_x, long cb_y, long cb_z, long chunk, int _unroll_p3);
void wave(float *  *  U_0_1_out, float *  U_0_m1, float *  U_0_0, float *  U_0_1, float fMin, float fDX, float fDT_DX_sq, long X_MAX, long Y_MAX, long Z_MAX);
void initialize_wave_parm(float *  U_0_m1, float *  U_0_0, float *  U_0_1, float fMin, float fDX, float fDT_DX_sq, long X_MAX, long Y_MAX, long Z_MAX, long cb_x, long cb_y, long cb_z, long chunk);
void initialize_wave(float *  U_0_m1, float *  U_0_0, float *  U_0_1, float fMin, float fDX, float fDT_DX_sq, long X_MAX, long Y_MAX, long Z_MAX);
