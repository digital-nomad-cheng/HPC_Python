echo "Generate data set with different grid size first.."
# loop over grid size
sizes=(60 120 240 480 960)
for size in "${sizes[@]}"
do 
  echo "grid size: $size"
  python3 generate-channel-grid.py --grid_size $size
  # loop over solvers
  solvers=("direct" "cg", "cgs", "gmres", "qmr", "lgmres")
  for solver in "${solvers[@]}"
  do
    echo "benchmarking solver: $solver"
    sed -e "s/solver_pp='direct'/solver_pp='$solver'/g" -e "s/solver_vel='direct'/solver_vel='$solver'/g" -e "s/solver_turb='direct'/solver_turb='$solver'/g" cupy_exec-pyCALC-RANS.py > cupy_test_$solver.py  
    sed -e "s/solver_pp='direct'/solver_pp='$solver'/g" -e "s/solver_vel='direct'/solver_vel='$solver'/g" -e "s/solver_turb='direct'/solver_turb='$solver'/g" numpy_exec-pyCALC-RANS.py > numpy_test_$solver.py  
    echo "run on cupy..."
    python3 -m cProfile -s tottime cupy_test_$solver.py | head -n 20 > cupy_${solver}_grid_size_${size}_out.txt
    echo "run on numpy..."
    python3 -m cProfile -s tottime numpy_test_$solver.py | head -n 20 > numpy_${solver}_grid_size_${size}_out.txt
  done
done


