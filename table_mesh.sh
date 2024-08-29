#!/bin/bash

# Variables
mesh_files=("project/mesh_test/hemisphere-quad-tri.obj" "project/mesh_test/hemisphere-cut.obj" "project/mesh_test/camelhead-cut.obj")
iterations=(20)
methods=("ARAP" "SYMMETRIC-DIRICHLET" "EXPONENTIAL-SYMMETRIC-DIRICHLET" "HENCKY-STRAIN" "AMIPS" "CONFORMAL-AMIPS-2D" "UNTANGLE-2D")

# Create the build directory if it doesn't exist
mkdir -p build

# Run cmake and make
cmake -B build -DCMAKE_BUILD_TYPE=Release && cd build && make -j
if [ $? -ne 0 ]; then
    echo "Build failed"
    exit 1
fi

for mesh in "${mesh_files[@]}"; do
    for iter in "${iterations[@]}"; do
        for method in "${methods[@]}"; do
            echo "Running project with mesh: $mesh, iterations: $iter, method: $method"
            if [ "$method" == "UNTANGLE-2D" ]; then
                ./project/SLIM_intern "$mesh" 1 max_iterations="$iter" "$method" epsilon=1e-1
            else
                ./project/SLIM_intern "$mesh" 1 max_iterations="$iter" "$method"
            fi
            if [ $? -ne 0 ]; then
                echo "Run failed with mesh: $mesh, iterations: $iter, method: $method"
            fi
        done
    done
done

echo "All tasks completed"