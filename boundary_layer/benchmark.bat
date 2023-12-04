@echo off
set "solvers=cgs"
set "gridsizes=1 10"

mkdir "bench_results" 2>nul

for %%g in (%gridsizes%) do (
    echo Benchmarking grid scalefactor: %%g
    python generate-bound-layer-grid.py %%g
    for %%s in (%solvers%) do (
        echo Benchmarking solver: %%s
        (
            python -m cProfile -s tottime cupy_exec-pyCALC-RANS.py %%s %%s %%s
        ) > "bench_results\cupy_%%s_grid_%%g_out.txt" 2>&1
    )
)