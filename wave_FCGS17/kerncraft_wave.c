float a[M][N][P];
float b[M][N][P];
float c[M][N][P];
float c1, c2, c3;

for(int i=2; i<M-2; ikerncraft -h
	++)
	for(int j=2; j<N-2; j++)
		for(int k=2; k<P-2; k++)
			c[i][j][k] = c1 * a[i][j][k] - b[i][j][k]
						+c2 * (a[i+1][j][k] + a[i-1][j][k] + a[i][j+1][k]
								+ a[i][j-1][k] + a[i][j][k+1] + a[i][j][k-1])
						+c3 * (a[i+2][j][k] + a[i-2][j][k] + a[i][j+2][k] 
								+ a[i][j-2][k] + a[i][j][k+2] + a[i][j][k-2]);

kerncraft -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml wave.c -p LC -p ECM -p Roofline -D M 200 -D N 200 -D P 200
                                   kerncraft                                    
wave.c                                  -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml
-D N 200 -D P 200 -D M 200
------------------------------------- ECM --------------------------------------
IACA analysis failed: pointer_increment could not be detected automatically. Use --pointer-increment to set manually to byte offset of store pointer address between consecutive assembly block iterations.
(myenv) ihpc07@lima1:~/stempel/stempel$ kerncraft -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml wave.c -p LC -p ECMData -p Roofline -D M 200 -D N 200 -D P 200
                                   kerncraft                                    
wave.c                                  -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml
-D P 200 -D M 200 -D N 200
----------------------------------- ECMData ------------------------------------
L1-L2 = 8.00 cy/CL
L2-L3 = 16.00 cy/CL
L3-MEM = 34.43 cy/CL

                                   kerncraft                                    
wave.c                                  -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml
-D P 200 -D M 200 -D N 200
-------------------------------------- LC --------------------------------------
1D layer condition:
L1: unconditionally fulfilled
L2: unconditionally fulfilled
L3: unconditionally fulfilled
2D layer condition:
L1: P <= 8192/11
L2: P <= 65536/11
L3: P <= 3145728/11
3D layer condition:
L1: 20*N*P + 8*P*(N - 2) + 16*P <= 32768
L2: 20*N*P + 8*P*(N - 2) + 16*P <= 262144
L3: 20*N*P + 8*P*(N - 2) + 16*P <= 12582912

                                   kerncraft                                    
wave.c                                  -m examples/machine-files/Intel_Xeon_CPU_X5650_2.67GHz_mod.yml
-D P 200 -D M 200 -D N 200
----------------------------------- Roofline -----------------------------------
Cache or mem bound with 1 core(s)
12.78 GFLOP/s due to MEM transfer bottleneck (bw with from copy benchmark)
Arithmetic Intensity: 1.00 FLOP/B


kerncraft -m examples/machine-files/AMD_Opteron_Processor_6274.yml wave.c -p LC -p ECMData -p Roofline -D M 200 -D N 200 -D P 200
                                   kerncraft                                    
wave.c                                  -m examples/machine-files/AMD_Opteron_Processor_6274.yml
-D N 200 -D M 200 -D P 200
-------------------------------------- LC --------------------------------------
Traceback (most recent call last):
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/bin/kerncraft", line 11, in <module>
    load_entry_point('kerncraft==0.6.0', 'console_scripts', 'kerncraft')()
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/lib/python3.5/site-packages/kerncraft/kerncraft.py", line 284, in main
    run(parser, args)
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/lib/python3.5/site-packages/kerncraft/kerncraft.py", line 250, in run
    model.analyze()
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/lib/python3.5/site-packages/kerncraft/models/layer_condition.py", line 225, in analyze
    self.results = self.calculate_cache_access()
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/lib/python3.5/site-packages/kerncraft/models/layer_condition.py", line 156, in calculate_cache_access
    csim = self.machine.get_cachesim(self._args.cores)
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/lib/python3.5/site-packages/kerncraft/machinemodel.py", line 66, in get_cachesim
    cs, caches, mem = cachesim.CacheSimulator.from_dict(cache_dict)
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/lib/python3.5/site-packages/cachesim/cache.py", line 63, in from_dict
    name=name, **{k:v for k,v in conf.items() if k not in ['store_to', 'load_from']})
  File "/home/hpc/ihpc/ihpc07/miniconda2/envs/myenv/lib/python3.5/site-packages/cachesim/cache.py", line 253, in __init__
    assert is_power2(ways), "ways needs to be a power of 2"
AssertionError: ways needs to be a power of 2
