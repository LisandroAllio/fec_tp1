class GaloisField:
    def __init__(self, m: int, pol: int):
        self.m: int = m
        self.pol: int = pol

    # ------------------- Validators ---------------------
    def validate(self, coef: int) -> bool:
        return coef < 0 or coef > (1 << self.m) - 1

    # ------------------- Operators ---------------------
    def suma(self, a: int, b: int) -> int:
        """La suma es sobre los coeficiente en mod 2. Analogo a una XOR"""
        return a ^ b

    def producto(self, a: int, b: int) -> int:
        """
        "Multiplicar dos elementos equivale a multiplicar sus polinomios y
        luego reducir el resultado módulo p(X). Esto es analogo a sumar (XOR)
        por p(X) en caso que este presente x^m.

        Con los binarios, por cada 1 de b (de menor a mayor) sumamos a. Luego
        multiplicamos a por x para pasar al siguiente bit y en caso de que el resultado
        contenga a x^m, reemplazamos
        """

        polFull = self.pol | (1 << self.m)
        r = 0

        while b > 0:
            if b & 1:
                r ^= a
            if a & (1 << self.m - 1):
                a = (a << 1) ^ polFull
            else:
                a = a << 1
            a &= (1 << self.m) - 1
            b >>= 1

        return r

    def inverso(self, a: int) -> int:
        """El inverso de a es a^-1, y a^-1 = a^(2^m - 2)"""
        if a == 0:
            raise ZeroDivisionError("El elemento 0 no tiene inverso multiplicativo")
        return self.potencia(a, (1 << self.m) - 2)  # a^(2^m - 2)

    def division(self, a: int, b: int) -> int:
        if b == 0:
            raise ZeroDivisionError("No se puede dividir por el elemento 0")
        return self.producto(a, self.inverso(b))

    def potencia(self, A: int, n: int) -> int:
        if n < 0:
            raise ValueError("El exponente debe ser no negativo")
        r = 1
        base = A
        while n > 0:
            if n & 1:
                r = self.producto(r, base)
            base = self.producto(base, base)
            n >>= 1
        return r

    def __eq__(self, otro):
        return self.m == otro.m and self.pol == otro.pol


class GFPoly:
    def __init__(self, gf: GaloisField, coefs: list):
        self._validate_pol(coefs, gf)
        self.gf: GaloisField = gf
        self.coefs: list = self._prune_zeros(coefs)

    @classmethod
    def contruir_con_raices(cls, gf: "GaloisField", raices: list) -> "GFPoly":
        f = cls(gf, [1])
        for r in raices:
            factor = cls(gf, [1, r])
            f = f * factor
        return f

    # ------------------- Validators ---------------------
    @staticmethod
    def _validate_pol(coefs: list, gf: GaloisField):
        for c in coefs:
            if gf.validate(c):
                raise ValueError(f"El coeficiente {c} no es un elemento del campo")

    @staticmethod
    def _prune_zeros(coefs: list) -> list:
        prune_coefs = list(coefs)
        while len(prune_coefs) > 1 and prune_coefs[0] == 0:
            prune_coefs.pop(0)
        return prune_coefs

    # ------------------- Operators ---------------------
    def __add__(self, otro: "GFPoly") -> "GFPoly":
        if self.gf != otro.gf:
            raise ValueError(
                "Los coeficientes de ambos polinomios deben pertenecer al mismo campo de Galois"
            )
        max_len = max(len(self.coefs), len(otro.coefs))
        r_coefs_rev = []
        for i in range(1, max_len + 1):
            c1 = self.coefs[-i] if i <= len(self.coefs) else 0
            c2 = otro.coefs[-i] if i <= len(otro.coefs) else 0
            r_coefs_rev.append(self.gf.suma(c1, c2))
        r_coefs = r_coefs_rev[::-1]
        return GFPoly(self.gf, r_coefs)

    def __mul__(self, otro: "GFPoly") -> "GFPoly":
        if self.gf != otro.gf:
            raise ValueError(
                "Los coeficientes de ambos polinomios deben pertenecer al mismo campo de Galois"
            )
        n = len(self.coefs) - 1
        k = len(otro.coefs) - 1
        r = [0] * (
            n + k + 1
        )  # la cantidad de coeficientes del resultado es (n + k + 1)
        for i in range(len(self.coefs)):
            for j in range(len(otro.coefs)):
                r[i + j] = self.gf.suma(
                    r[i + j], self.gf.producto(self.coefs[i], otro.coefs[j])
                )
        return GFPoly(self.gf, r)

    def __floordiv__(self, otro: "GFPoly") -> list["GFPoly"]:
        if self.gf != otro.gf:
            raise ValueError(
                "Los coeficientes de ambos polinomios deben pertenecer al mismo campo de Galois"
            )
        if len(otro.coefs) == 1 and otro.coefs[0] == 0:
            raise ZeroDivisionError("No se puede dividir por un polinomio nulo")

        n = len(self.coefs) - 1
        k = len(otro.coefs) - 1

        cociente_x: list = [0] * (n - k + 1)
        g = len(cociente_x) - 1
        resto_x: list = list(self.coefs)

        while (len(resto_x) - 1) >= k:
            r = len(resto_x) - 1
            c_r = resto_x[0]
            c_g = otro.coefs[0]

            t = self.gf.division(c_r, c_g)
            m = r - k

            termino = [self.gf.producto(coef, t) for coef in otro.coefs]
            termino.extend([0] * m)

            resto_x = self._prune_zeros(
                [self.gf.suma(resto_x[i], termino[i]) for i in range(len(resto_x))]
            )
            cociente_x[g - m] = t

        return [GFPoly(self.gf, cociente_x), GFPoly(self.gf, resto_x)]

    def escalado(self, factor: int) -> "GFPoly":
        self.gf.validate(factor)
        r = [self.gf.producto(coef, factor) for coef in self.coefs]
        return GFPoly(self.gf, r)

    # Alfortirmo de Horner:
    # a_n·X^n + a_(n-1)·X^(n-1) + ... + a_1·X + a_0 = ( ... ( ( a_n·X + a_(n-1) )·X + a_(n-2) )·X + ... + a_1 )·X + a_0
    #
    def evaluar(self, valor: int) -> int:
        r = self.coefs[0]
        for c in self.coefs[1:]:
            r = self.gf.suma(self.gf.producto(r, valor), c)
        return r
