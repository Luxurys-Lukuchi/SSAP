import math


def microstrip_corner(w_um: float, h_um: float, er: float):

    w = w_um * 1e-6          # перевод в метры
    h = h_um * 1e-6
    wh = w / h              

    # Формула 11.84
    C = w * ((10.35 * er + 2.5) * wh + (2.6 * er + 5.64))       # пФ

    # Формула 11.85
    L = 220.0 * h * (1.0 - 1.35 * math.exp(-0.18 * wh ** 1.39))  # нГн

    return C, L


if __name__ == "__main__":
    W, H, ER = 50.0, 100.0, 12.9   

    C, L = microstrip_corner(W, H, ER)

    print(f"Исходные данные: W = {W} мкм, h = {H} мкм, Er = {ER}")
    print(f"W/h = {W / H:.3f}")
    print(f"C = {C:.6f} пФ  = {C * 1e3:.3f} фФ")
    print(f"L = {L:.6f} нГн = {L * 1e3:.3f} пГн")


import math

ZF0 = 376.73  


def microstrip_line(w_um, h_um, er):
   
    u = w_um / h_um
    a = (1.0 + math.log((u**4 + (u / 52.0)**2) / (u**4 + 0.432)) / 49.0
         + math.log(1.0 + (u / 18.1)**3) / 18.7)
    b = 0.564 * ((er - 9.0) / (er + 3.0)) ** 0.053
    eeff = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)
    F = 6.0 + (2 * math.pi - 6.0) * math.exp(-(30.666 / u) ** 0.7528)
    z0 = ZF0 / (2 * math.pi * math.sqrt(eeff)) * math.log(F / u + math.sqrt(1 + (2 / u)**2))
    return z0, eeff


def microstrip_cross(w1_um, w2_um, h_um, er, improve=True):
    """Крест: C1..C4 [Ф], L1..L3 [Гн]. W1, W2, h — в мкм."""
    w1, w2, h = w1_um * 1e-6, w2_um * 1e-6, h_um * 1e-6
    u1, u2 = w1 / h, w2 / h

    # 11.226, 11.227 (база, er = 9.9)
    X = math.log10(u1) * (86.6 * u2 - 30.9 * math.sqrt(u2) + 367) + u2**3 + 74 * u2 + 130
    C = 1e-12 * w1 * (0.25 * X * (1 / u1)**(1/3) - 60 + 1 / (2 * u2) - 0.375 * u1 * (1 - u2))

    # 11.228, 11.229;  L2 — обмен индексов W1 <-> W2
    Y = 165.6 * u2 + 31.2 * math.sqrt(u2) - 11.8 * u2**2
    L1 = 1e-9 * h * (Y * u1 - 32 * u2 + 3) * (1 / u1)**1.5
    Y2 = 165.6 * u1 + 31.2 * math.sqrt(u1) - 11.8 * u1**2
    L2 = 1e-9 * h * (Y2 * u2 - 32 * u1 + 3) * (1 / u2)**1.5

    # 11.230
    L3 = 1e-9 * h * (5 * u2 * math.cos(math.pi / 2 * (1.5 - u1)) - (1 + 7 / u1) / u2 - 337.5)

    C1 = C2 = C3 = C4 = C
    if improve:
        L3 *= 0.8                                   

        if er != 9.9:                               
            za, ea = microstrip_line(w1_um, h_um, 9.9)
            zb, eb = microstrip_line(w1_um, h_um, er)
            k1 = za / zb * math.sqrt(eb / ea)
            za, ea = microstrip_line(w2_um, h_um, 9.9)
            zb, eb = microstrip_line(w2_um, h_um, er)
            k2 = za / zb * math.sqrt(eb / ea)
            C1 = C3 = C * k1                        
            C2 = C4 = C * k2                        

    return {"C1": C1, "C2": C2, "C3": C3, "C4": C4, "L1": L1, "L2": L2, "L3": L3}


if __name__ == "__main__":
    W1, W2, H, ER = 50.0, 50.0, 100.0, 12.9
    r = microstrip_cross(W1, W2, H, ER)
    print(f"W1 = {W1} мкм, W2 = {W2} мкм, h = {H} мкм, Er = {ER}")
    for k in ("C1", "C2", "C3", "C4"):
        print(f"{k} = {r[k]*1e12:.6f} пФ = {r[k]*1e15:.3f} фФ")
    for k in ("L1", "L2", "L3"):
        print(f"{k} = {r[k]*1e9:.6f} нГн = {r[k]*1e12:.3f} пГн")
