def sign_changed_here(X):
	import numpy as np
	sign_change_index = []
	if np.shape(X) == (len(X),): 
		Xsign = np.sign(X[len(X):])
		signchange = ((np.roll(Xsign,1)-Xsign) != 0).astype(int)
		signchange[0] = 0
		if 1 in list(signchange): 
			while 1 in list(signchange):
				sign_change_index.append(list(signchange).index(1))
				signchange[list(signchange).index(1)]=0
		else:
			sign_change_index.append('No Sign Change Found!')
	else:
		for i in range(np.shape(X)[0]):
			Xsign = np.sign(X[i,:])
			Xsign[0] = Xsign[1]
			signchange = ((np.roll(Xsign,1)-Xsign) != 0).astype(int)
			signchange[0] = 0
			if 1 in list(signchange):
				temp = []
				while 1 in list(signchange): 
					temp.append(list(signchange).index(1))
					signchange[list(signchange).index(1)] = 0 
				sign_change_index.append(temp[-1])
			else:
				sign_change_index.append('No Sign Change Found!')
	return(sign_change_index)

