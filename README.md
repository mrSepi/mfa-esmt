A heuristic algorithm for the Euclidean Steiner Minimal Tree.
- `mfa-mst` doesn't require any external library
- `mfa-mst-delaunay` uses [Triangle](https://people.math.sc.edu/Burkardt/c_src/triangle/triangle.html) to generate Delaunay triangulations to speed up the computation
- `mfa-lp` uses [GLPK](https://www.gnu.org/software/glpk/) to solve linear programming problems

The version `mfa-lp` has exponential time complexity, while `mfa-mst` cubic.
Given as input a list of 2-dimensional points in the form
```
<number of points>
 <x> <y>
...
```
returns the coordinates of the Steiner points.
