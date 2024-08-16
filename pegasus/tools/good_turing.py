import numpy as np
from scipy.stats import linregress


class SimpleGoodTuring():
    def __init__(self, counts):
        r, Nr = np.unique(counts, return_counts=True)
        if r[0] == 0:
            n0 = Nr[0]
            Nr = Nr[1:]
            r = r[1:]
        else:
            n0 = 0
        self._counts = counts
        self._r = r
        self._Nr = Nr
        self._n0 = n0

        self.smooth_fit()
        self.estimate()

    def smooth_fit(self):
        z = np.zeros(self._r.size)
        for i in range(z.size):
            if i == 0:
                z[i] = 2 * self._Nr[i] / self._r[i + 1]
            elif i == z.size - 1:
                z[i] = self._Nr[i] / (self._r[i] - self._r[i - 1])
            else:
                z[i] = 2 * self._Nr[i] / (self._r[i + 1] - self._r[i - 1])

        # Perform linear regression on log-log scale
        x = np.log(self._r)
        y = np.log(z)
        slope, intercept, _, _, _ = linregress(x, y)
        self._slope = slope
        self._intercept = intercept

    def smoothed_Nr(self, r):
        return np.exp(self._intercept + self._slope * np.log(r))

    def std_r(self, r):
        i = np.where(self._r==r)[0][0]
        return np.sqrt((r + 1)**2 * (self._Nr[i + 1] / self._Nr[i]**2) * (1 + self._Nr[i + 1] / self._Nr[i]))

    def estimate(self):
        # Calculate p0
        N = np.sum(self._r * self._Nr)
        p0 = self._Nr[0] / N if 1 in self._r else 0

        # Calculate Good-Turing estimates
        r_star = np.zeros(self._r.size)
        for i, cur_r in enumerate(self._r):
            y_r = (cur_r + 1) * self.smoothed_Nr(cur_r + 1) / self.smoothed_Nr(cur_r)
            if (cur_r + 1) in self._r:
                x_r = (cur_r + 1) * self._Nr[i + 1] / self._Nr[i]
                r_star[i] = y_r if np.abs(x_r - y_r) <= 1.96 * self.std_r(cur_r) else x_r
            else:
                r_star[i] = y_r

        N_estimated = np.sum(r_star * self._Nr)
        prob_estimated = np.zeros(r_star.size)
        for i, cur_r_star in enumerate(r_star):
            prob_estimated[i] = (1 - p0) * cur_r_star / N_estimated
        prob_unseen = p0 / self._n0 if self._n0 > 0 else 0

        prop_adj = np.zeros(self._counts.size)
        for i in range(prop_adj.size):
            cur_r = self._counts[i]
            if cur_r > 0:
                j = np.where(self._r==cur_r)[0][0]
                prop_adj[i] = prob_estimated[j]
            else:
                prop_adj[i] = prob_unseen

        self._prop_adj = prop_adj
        self._p0 = p0

    def get_proportions(self):
        return self._prop_adj
