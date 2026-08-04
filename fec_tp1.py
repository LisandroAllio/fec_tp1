class GaloiField:
    def __init_(self, m, pol):
        self.m: int = m
        self.pol: int = pol

    def sum(self, a: int, b: int) -> int:
        return a ^ b

    def producto(self, a: int, b: int) -> int:
        polFull = self.pol | (1 << self.m)
        r = 0

        while b > 0:
            if b & 1:
                r += a
            if a & (1 << self.m-1):
                a = (a << 1) ^ polFull
            else:
                a = (a << 1)
            a &= (1<< self.m)-1
            b >>= 1

        return r

    def inverso(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("El elemento 0 no tiene inverso multiplicativo")
        return self.potencia(a, (1 << self.m) - 2)   # a^(2^m - 2)

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
