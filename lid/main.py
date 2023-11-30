from data import generate_input_data
from cupy_exec_pyCALC_RANS import test_cupy_solver
import numpy as np 

if __name__ == "__main__":
    generate_input_data(120)
    test_cupy_solver()

