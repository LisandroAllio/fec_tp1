# fec_tp1

**Trabajo Práctico 1** del curso de FEC 2026 (Fundación Fulgor) — **Álgebra en FEC.**

## Consigna

El TP consiste en implementar en Python estructuras algebraicas sobre campos de Galois `GF(2^m)`:

1. **Campos de Galois**: una clase que representa `GF(2^m)`, construida a partir del orden `m`
   y el polinomio primitivo `P(x)` (como entero de `m` bits, sin el término `x^m`). Debe soportar
   suma, producto, inverso multiplicativo, división y potencia entre elementos del campo
   (representados como enteros en `[0, 2^m - 1]`).
2. **Polinomios sobre un campo de Galois** (`GFPoly`): una clase que representa polinomios con
   coeficientes en `GF(2^m)`, soportando suma, producto, división entera (cociente y resto),
   escalado, evaluación y construcción a partir de raíces.

El enunciado completo está en [`TP1_FEC.pdf`](TP1_FEC.pdf).
