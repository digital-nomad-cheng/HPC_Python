1. replace numpy with cupy 
2. replace np.insert with cp.append
3. replace np.matrix.flatten(a) with a.flatten()


# How to run
1. profile program sort function by total time it consumes: 
	`python3 -m cProfile -s tottime cupy_exec-pyCALC-RANS.py > profiling_results.txt`
    `python -m cProfile -s tottime your_script.py | head -n 20`
