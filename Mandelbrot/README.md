This projects provides two implementations of Mandelbrot set using GPU.

In Mandelbrot1 every GPU block is used to compute 1 row of a pixel matrix.
Mandelbrot2 uses better implementation, where every block is used to compute values in a 8x8 submatrix.
