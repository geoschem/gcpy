#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Some Quantum Mechanics, filling an atomic orbital
=================================================

Considering an atomic single orbital and how to fill it by use of the
chemical potential.
"""
# Code source: Óscar Nájera
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np
mu = np.linspace(0, 3, 800)
for b in [10, 20, 30]:
    n = 2 * (np.exp(b * (mu - 1)) + np.exp(b * (2 * mu - 3))) / \
        (1 + np.exp(b * (mu - 1)) * (2 + np.exp(b * (mu - 2))))
    plt.plot(mu, n, label=r"$\beta={0}$".format(b))
plt.xlabel(r'$\mu$ ($\epsilon=1$, $U=1$)')
plt.ylabel(r'$\langle N \rangle=\langle n_\uparrow \rangle+\langle n_\downarrow\rangle$')
plt.legend(loc=0)
plt.show()
