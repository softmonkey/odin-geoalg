# Geometric Algebra for Odin

An implementation of [Geometric Algebra](https://bivector.net/) for in the [Odin](https://odin-lang.org/) language. Other algebras will come as needed - issues welcome.

## Supported Algebras
* [3D Projective Geometric Algebra](https://bivector.net/tools.html?p=3&q=0&r=1), a.k.a. PGA3D
* [2D Projective Geometric Algebra](https://bivector.net/tools.html?p=2&q=0&r=1), a.k.a. PGA2D

Currently available for f32. Plan to add support for f16 and f64 later - submit an issue if needed.

## Notes
* Formating of multi-vectors uses unicode for basis vectors and colors for additional clarity (on dark backgrounds at least). For example:
  ```powershell
  point on plane: [0.200e₀₂₁, 1.400e₀₃₂, 1.000e₁₂₃]
  ```
  These features can be disabled if prefered then we prefer to avoid confusion with the exponential notation. Without Unicode, formating a multi-vector will result in:  
  ```
  [1.000e0, 2.000e1, -1.500e2, -3.000e0e1]
  ```
  The use of unicode is just for pretty formating, the basis vector constants for each algebra are addressed as regular characters.

## Examples
Coming soon!
