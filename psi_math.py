import numpy as np

# this file contains "full" expressions 
# for the instantaneous bottleneck model,
# without any mathematical simplification.
# these functions are used in plotting


# TOPOLOGY PROBABILITIES

def p_alpha(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB

    ans = (
        np.exp(-t / nB - t / nA - s / nA) / 9
        + (1 - np.exp(-t / nB)) * np.exp(-t / nA - s / nA) / 3
        + np.exp(-t / nB) * (1 - np.exp(-t / nA - s / nA)) / 3
        + (1 - np.exp(-t / nB)) * (1 - np.exp(-t / nA - s / nA))
    )

    return ans


def p_beta_B(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB
    # cherry in B
    ans = (2 / 3) * np.exp(-t / nA - s / nA) * (1 - np.exp(-t / nB))

    return ans


def p_beta_C(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB
    # cherry in C
    ans = (1 / 9) * np.exp(-t / nA - s / nA) * np.exp(-t / nB)

    return ans


def p_gamma_A(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB

    # cherry in A
    ans = (2 / 3) * (1 - np.exp(-t / nA)) * np.exp(-t / nB)

    return ans


def p_gamma_C(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB

    # cherry in C
    ans = (1 / 9) * np.exp(-t / nA - s / nA) * np.exp(-t / nB)

    return ans


def p_gamma_bot(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB

    # cherry in bottleneck
    ans = (2 / 3) * np.exp(-t / nA - t / nB) * (1 - np.exp(-s / nA))

    return ans


def p_delta(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB

    ans = 2 * np.exp(-t / nA - t / nB - s / nA) / 9

    return ans


def p_eps(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB

    ans = 2 * np.exp(-t / nA - t / nB - s / nA) / 9

    return ans


def p_zeta(s, t, nA, nB):
    nA = 2 * nA
    nB = 2 * nB

    ans = 2 * np.exp(-t / nA - t / nB - s / nA) / 9

    return ans


# expected total lengths


def l_beta_C(s, t, nA, nB, nC):
    # cherry in c
    return 4 * t + 22 * nC / 3


def l_beta_B(s, t, nA, nB, nC):
    return 2 * nB + 3 * t + 6 * nC - t / (np.exp(t / 2 / nB) - 1)


def l_gamma_C(s, t, nA, nB, nC):
    return 4 * t + 22 * nC / 3


def l_gamma_A(s, t, nA, nB, nC):
    return 2 * nA + 3 * t + 6 * nC - t / (np.exp(t / 2 / nA) - 1)


def l_gamma_bot(s, t, nA, nB, nC):
    return 4 * t + 6 * nC


def l_delta(t, nC):
    return 4 * t + 22 * nC / 3


def l_eps(t, nC):
    return 4 * t + 22 * nC / 3


def l_zeta(t, nC):
    return 4 * t + 22 * nC / 3


# expected branch lengths


def b_beta(nC):
    return 2 * nC


def b_gamma(nC):
    return 2 * nC


def b_epsilon_1(nC):
    return 2 * nC / 3


def b_epsilon_2(nC):
    return 2 * nC


def b_zeta_1(nC):
    return 2 * nC / 3


def b_zeta_2(nC):
    return 2 * nC


def b_delta(nC):
    return 4 * nC + 2 * nC / 3


# now we need shared branch lengths p11, p12, p21


def p11(s, t, nA, nB, nC):
    ans = p_delta(s, t, nA, nB) * b_delta(nC)
    ans += p_eps(s, t, nA, nB) * b_epsilon_1(nC)
    ans += p_zeta(s, t, nA, nB) * b_zeta_1(nC)
    return ans


def p12(s, t, nA, nB, nC):
    ans = p_beta_B(s, t, nA, nB) * b_beta(nC) + p_beta_C(s, t, nA, nB) * b_beta(nC)
    ans += p_eps(s, t, nA, nB) * b_epsilon_2(nC)
    return ans


def p21(s, t, nA, nB, nC):
    ans = p_gamma_A(s, t, nA, nB) * b_gamma(nC) + p_gamma_C(s, t, nA, nB) * b_gamma(nC)
    ans += p_gamma_bot(s, t, nA, nB) * b_gamma(nC)
    ans += p_zeta(s, t, nA, nB) * b_zeta_2(nC)
    return ans


# psi

def epsi(s, t, nA, nB, nC):
    return (p21(s, t, nA, nB, nC) - p12(s, t, nA, nB, nC)) / (
        p11(s, t, nA, nB, nC) + p21(s, t, nA, nB, nC) + p12(s, t, nA, nB, nC)
    )


def epsisq(s, t, nA, nB, nC):
    return (p21(s, t, nA, nB, nC) + p12(s, t, nA, nB, nC)) / (
        p11(s, t, nA, nB, nC) + p21(s, t, nA, nB, nC) + p12(s, t, nA, nB, nC)
    )


def varpsi(s, t, nA, nB, nC):
    return epsisq(s, t, nA, nB, nC) - epsi(s, t, nA, nB, nC) ** 2
