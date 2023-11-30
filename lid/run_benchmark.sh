solvers=("direct" "cg")
for solver in "${solvers[@]}"
do
  echo "benchmarking solver: $solver" 
  sed -e "s/solver_pp='direct'/solver_pp='$solver'/g" -e "s/solver_vel='direct'/solver_vel='$solver'/g" -e "s/solver_turb='direct'/solver_turb='$solver'/g" cupy_exec-pyCALC-RANS.py > test_$solver.py  
  python3 -m cProfile -s tottime test_$solver.py | head -n 20 > cupy_${solver}_grid_120_out.txt
done
# sed -i 's/old_string/new_string/g' filename.txt

