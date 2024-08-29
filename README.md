# Compile and run:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release && cd build && make -j && project/SLIM_intern name=$(Name) weights=$(Weights) max_iterations=$(Max Iterations) energy=$(Energy) epsilon=$(UNTANGLE 2D)
```

### Detailed Parameter Choices

1. **Name**: Path to your mesh file
   - Default: `project/mesh_test/hemisphere.obj`

2. **Weights**: Type of weights for Tutte's embedding
   - Default: `1` (uniform weights)
   - Options:
     - `1`: Uniform weights
     - `2`: Cotangent weights
     - `3`: Random weights
     - `4`: Swap points
     - `5`: Texture coordinates, only extract the texture coordinates of the mesh

3. **Energy**: Distortion measure used in the optimization
   - Default: `"ARAP"`
   - Options:
     - `"ARAP"`
     - `"SYMMETRIC_DIRICHLET"`
     - `"EXPONENTIAL_SYMMETRIC_DIRICHLET"`
     - `"HENCKY_STRAIN"`
     - `"AMIPS"`
     - `"CONFORMAL_AMIPS_2D"`
     - `"UNTANGLE_2D"`

4. **Max Iterations**: Maximum number of iterations for the optimization
   - Default: `100`
   - Different value in the format `max_iterations=VALUE`.

5. **UNTANGLE 2D**: Specify the value of the initial epsilon.
   - Default: `1e-1`
   - Different value in the format `epsilon=VALUE`.

### Example Usage

You can provide the arguments directly. For instance:

```sh
make -j && project/SLIM_intern name=mesh_test/hemisphere.obj weights=1 max_iterations=20 energy=UNTANGLE_2D epsilon=0.5
```

This command uses:
- `mesh_test/hemisphere.obj` as the mesh file.
- `1` for uniform weights.
- `20` as the maximum number of iterations.
- `UNTANGLE_2D` as the energy.
- `0.5` as the initial epsilon value for untangling in 2D.

### References

DOI:10.1145/3450626.3459847

DOI:10.1145/2983621
