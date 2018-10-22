def estimate_probs_old(arr, pvec, n_iter = 200):
	probs = np.zeros(pvec.size + 1)
	z = np.zeros(pvec.size + 1)
	noise = pvec.size
	probs[:] = 0.1 / pvec.size
	probs[-1] = 0.9
	for i in range(n_iter):
		# E step
		z[:] = 0.0
		for j in range(pvec.size):
			if arr[j] > 0:
				p = probs[j] / (probs[noise] * pvec[j] + probs[j])
				z[j] += arr[j] * p
				z[noise] += arr[j] * (1.0 - p)
		# M step
		probs = z / z.sum()
	return probs

def remove_ambient(arr, prob_noise, pvec, n_iter = 50):
	idx = arr > 0.0
	size = sum(idx)

	newarr = arr[idx]
	probs = np.zeros(size)
	probs[:] = (1.0 - prob_noise) / size
	nprobs = pvec[idx] / pvec[idx].sum()
	print(probs)
	# EM algorithm
	z = np.zeros(size)
	for i in range(n_iter):
		# E step
		z[:] = 0.0
		for j in range(size):
			p = probs[j] / (prob_noise * nprobs[j] + probs[j])
			z[j] += newarr[j] * p
		# M step
		old_probs = probs
		probs = (1.0 - prob_noise) * z / z.sum()
		print("Diff = {}".format(np.linalg.norm(probs - old_probs, ord = 1)))
		print(probs)
		print(sum(probs > 0))

	results = np.zeros(pvec.size)
	results[idx] = probs

	return results

