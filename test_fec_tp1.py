import unittest

from fec_tp1 import GaloisField, GFPoly


class TestGaloisField(unittest.TestCase):
    def setUp(self):
        # GF(2^3) con p(x) = x^3 + x + 1 (pol = 0b011, sin el bit x^3 implicito)
        self.gf = GaloisField(3, 0b011)

    def test_suma_es_xor(self):
        self.assertEqual(self.gf.suma(3, 5), 3 ^ 5)
        self.assertEqual(self.gf.suma(0, 7), 7)
        self.assertEqual(self.gf.suma(6, 6), 0)

    def test_producto_dentro_del_campo(self):
        for a in range(8):
            for b in range(8):
                r = self.gf.producto(a, b)
                self.assertGreaterEqual(r, 0)
                self.assertLessEqual(r, 7)

    def test_producto_valores_conocidos(self):
        # (x+1)(x^2+1) = x^3+x^2+x+1 -> x^3 = x+1 -> reduce a x^2 = 4
        self.assertEqual(self.gf.producto(3, 5), 4)
        self.assertEqual(self.gf.producto(0, 5), 0)
        self.assertEqual(self.gf.producto(1, 5), 5)

    def test_producto_es_conmutativo(self):
        for a in range(8):
            for b in range(8):
                self.assertEqual(self.gf.producto(a, b), self.gf.producto(b, a))

    def test_inverso_multiplicativo(self):
        for a in range(1, 8):
            inv = self.gf.inverso(a)
            self.assertEqual(self.gf.producto(a, inv), 1)

    def test_inverso_de_cero_lanza_excepcion(self):
        with self.assertRaises(ZeroDivisionError):
            self.gf.inverso(0)

    def test_division(self):
        for a in range(8):
            for b in range(1, 8):
                # (a / b) * b == a
                cociente = self.gf.division(a, b)
                self.assertEqual(self.gf.producto(cociente, b), a)

    def test_division_por_cero_lanza_excepcion(self):
        with self.assertRaises(ZeroDivisionError):
            self.gf.division(1, 0)

    def test_potencia(self):
        self.assertEqual(self.gf.potencia(5, 0), 1)
        self.assertEqual(self.gf.potencia(5, 1), 5)
        self.assertEqual(self.gf.potencia(5, 2), self.gf.producto(5, 5))
        self.assertEqual(self.gf.potencia(5, 3), self.gf.producto(self.gf.producto(5, 5), 5))

    def test_potencia_exponente_negativo_lanza_excepcion(self):
        with self.assertRaises(ValueError):
            self.gf.potencia(5, -1)

    def test_validate(self):
        # validate() devuelve True cuando el coeficiente esta FUERA de rango
        self.assertFalse(self.gf.validate(0))
        self.assertFalse(self.gf.validate(7))
        self.assertTrue(self.gf.validate(-1))
        self.assertTrue(self.gf.validate(8))

    def test_eq(self):
        otro = GaloisField(3, 0b011)
        distinto = GaloisField(4, 0b011)
        self.assertEqual(self.gf, otro)
        self.assertNotEqual(self.gf, distinto)


class TestGFPoly(unittest.TestCase):
    def setUp(self):
        self.gf = GaloisField(3, 0b011)

    def test_prune_zeros_en_construccion(self):
        p = GFPoly(self.gf, [0, 0, 3, 5])
        self.assertEqual(p.coefs, [3, 5])

    def test_prune_zeros_deja_al_menos_un_coeficiente(self):
        p = GFPoly(self.gf, [0, 0, 0])
        self.assertEqual(p.coefs, [0])

    def test_constructor_valida_coeficientes(self):
        with self.assertRaises(ValueError):
            GFPoly(self.gf, [1, 8])  # 8 esta fuera del campo GF(2^3)

    def test_suma(self):
        p1 = GFPoly(self.gf, [1, 2, 3])  # x^2 + 2x + 3
        p2 = GFPoly(self.gf, [4, 5])  # 4x + 5
        r = p1 + p2
        self.assertEqual(r.coefs, [1, self.gf.suma(2, 4), self.gf.suma(3, 5)])

    def test_suma_campos_distintos_lanza_excepcion(self):
        otro_gf = GaloisField(4, 0b011)
        p1 = GFPoly(self.gf, [1])
        p2 = GFPoly(otro_gf, [1])
        with self.assertRaises(ValueError):
            p1 + p2

    def test_mul_por_cero_da_cero(self):
        p1 = GFPoly(self.gf, [3, 5])
        cero = GFPoly(self.gf, [0])
        r = p1 * cero
        self.assertEqual(r.coefs, [0])

    def test_mul_por_uno_es_identidad(self):
        p1 = GFPoly(self.gf, [3, 5, 2])
        uno = GFPoly(self.gf, [1])
        r = p1 * uno
        self.assertEqual(r.coefs, p1.coefs)

    def test_mul_grado_resultante(self):
        p1 = GFPoly(self.gf, [1, 0])  # x
        p2 = GFPoly(self.gf, [1, 0])  # x
        r = p1 * p2  # x^2
        self.assertEqual(r.coefs, [1, 0, 0])

    def test_floordiv_cociente_y_resto(self):
        # (x^2 + 2x + 3) / (x + 1) -> verificar que d*q + r == a
        a = GFPoly(self.gf, [1, 2, 3])
        b = GFPoly(self.gf, [1, 1])
        cociente, resto = a // b
        reconstruido = (b * cociente) + resto
        self.assertEqual(reconstruido.coefs, a.coefs)

    def test_floordiv_por_polinomio_nulo_lanza_excepcion(self):
        a = GFPoly(self.gf, [1, 2, 3])
        cero = GFPoly(self.gf, [0])
        with self.assertRaises(ZeroDivisionError):
            a // cero

    def test_escalado(self):
        p = GFPoly(self.gf, [3, 5])
        factor = 2
        r = p.escalado(factor)
        esperado = [self.gf.producto(3, factor), self.gf.producto(5, factor)]
        self.assertEqual(r.coefs, esperado)

    def test_evaluar_constante(self):
        p = GFPoly(self.gf, [4])
        self.assertEqual(p.evaluar(0), 4)
        self.assertEqual(p.evaluar(6), 4)

    def test_evaluar_lineal(self):
        # p(x) = x + 3 -> p(2) = 2 + 3 (suma en GF, no aritmetica normal)
        p = GFPoly(self.gf, [1, 3])
        self.assertEqual(p.evaluar(2), self.gf.suma(2, 3))

    def test_construir_con_raices(self):
        raices = [2, 3]
        p = GFPoly.contruir_con_raices(self.gf, raices)
        # el polinomio construido debe anularse en cada raiz
        for r in raices:
            self.assertEqual(p.evaluar(r), 0)


if __name__ == "__main__":
    unittest.main()
