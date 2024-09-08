# Compile and run:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release &&
cd build &&
make -j &&
project/SLIM_intern name=project/mesh_test/hemisphere.obj weights=1 max_iterations=20 energy=UNTANGLE_2D epsilon=0.5
```

This command uses:
- `project/mesh_test/hemisphere.obj` as the mesh file.
- `1` for uniform weights.
- `20` as the maximum number of iterations.
- `UNTANGLE_2D` as the energy.
- `0.5` as the initial epsilon value for untangling in 2D.

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
   - Default: `"SYMMETRIC-DIRICHLET"`
   - Options:
     - `"ARAP"`
     - `"SYMMETRIC-DIRICHLET"`
     - `"EXPONENTIAL-SYMMETRIC-DIRICHLET"`
     - `"HENCKY-STRAIN"`
     - `"AMIPS"`
     - `"CONFORMAL-AMIPS-2D"`
     - `"UNTANGLE-2D"`

4. **Max Iterations**: Maximum number of iterations for the optimization
   - Default: `100`
   - Different value in the format `max_iterations=VALUE`.

5. **For `"UNTANGLE_2D"` Only**: Specify the value of the initial epsilon.
   - Default: `1e-1`
   - Different value in the format `epsilon=VALUE`.

### References

Vladimir Garanzha, Igor Kaporin, Liudmila Kudryavtseva, François Protais, Nicolas Ray, and Dmitry Sokolov. 2021. Foldover-free maps in 50 lines of code. ACM Trans. Graph. 40, 4, Article 102 (August 2021), 16 pages. https://doi.org/10.1145/3450626.3459847

Michael Rabinovich, Roi Poranne, Daniele Panozzo, and Olga Sorkine-Hornung. 2017. Scalable Locally Injective Mappings. ACM Trans. Graph. 36, 2, Article 16 (April 2017), 16 pages. https://doi.org/10.1145/2983621
