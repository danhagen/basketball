import ipdb
import numpy as np 
import matplotlib.pyplot as plt
import pickle

EndTime = 0.55
ChangeInTime = 0.0001
Time = np.arange(0,EndTime+ChangeInTime,ChangeInTime,float)

from GenerateTrajectories import Spline

Angle1Splines = pickle.load(open('SplineClassObjects.pkl','rb'))[0]
Angle2Splines = pickle.load(open('SplineClassObjects.pkl','rb'))[1]
Angle3Splines = pickle.load(open('SplineClassObjects.pkl','rb'))[2]

def plot_spline_results(Time,AngleSplines):
	"""
	For a 1D vector Time and an N-D list of arrays AngleSplines, this will plot all calculated splines on the
	same figure. Must close final figure in order to run next line of code.
	"""
	plt.figure()
	for i in range(0,AngleSplines.shape[0]):
		plt.plot(Time,AngleSplines[i].pp_func(Time))
	plt.show()
def normalized_muscle_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
	#R_tranpose Column 1
	r1DELTa = 19
	r1CB = 20
	r1DELTp = -8
	r1BIC = 15
	r1TRI = -15
	r1BRA = 0
	r1BRD = 0
	r1PRO = 0
	r1FCR = 0
	r1ECRB = 0
	r1ECRL = 0
	r1FCU = 0
	r1FDS = 0
	r1PL = 0
	r1ECU = 0
	r1EDM = 0
	r1EDC = 0
	r1APL = 0

	#R_tranpose Column 2
	r2DELTa = 0
	r2CB = 0
	r2DELTp = 0
	r2BIC = lambda Angle2: (24777*np.pi)/10000 + (1288265228720957*Angle2)/35184372088832 \
							- (2429*np.pi*Angle2)/125 + (68251*np.pi*(Angle2**2))/5000 \
							- (10427*np.pi*(Angle2**3))/5000 + (20571*(Angle2**2))/10000 \
							- (14043*(Angle2**3))/2500 + 84533/10000
	r2TRI = lambda Angle2: - (8759*(Angle2**3))/5000 + (93509*(Angle2**2))/10000 \
							- (88691*Angle2)/10000 - 863614486669217/35184372088832
	r2BRA = lambda Angle2: - (12667*(Angle2**3))/2000 + (30689*(Angle2**2))/1250 \
							- (4544779416463265*Angle2)/281474976710656 + 1139910323808397/70368744177664
	r2BRD = lambda Angle2: (28129*np.pi)/10000 - (23671*Angle2)/2000 - (57781*np.pi*Angle2)/10000 \
							+ (3629*np.pi*(Angle2**2))/1250 - (197*np.pi*(Angle2**3))/500 \
							+ (24636921970321*(Angle2**2))/549755813888 - (33739*(Angle2**3))/2500 + 38141/2500
	r2PRO = lambda Angle2: (3933*np.pi)/10000 - (10079*Angle2)/10000 - (13103*np.pi*Angle2)/1250 \
							+ (2597831076304493*np.pi*(Angle2**2))/70368744177664 + (2202*np.pi**2*Angle2)/625 \
							- (93111*np.pi*(Angle2**3))/2500 + (72987*np.pi*(Angle2**4))/5000 \
							- (20089*np.pi*(Angle2**5))/10000 - (4369*np.pi**2)/10000 \
							- (6847666938421497*(Angle2**2))/562949953421312 + (53151*(Angle2**3))/2500 \
							- (5503*(Angle2**4))/500 + (8763*(Angle2**5))/5000 \
							- (1466808324885735*np.pi**2*(Angle2**2))/140737488355328 + (51333*np.pi**2*(Angle2**3))/5000 \
							- (39919*np.pi**2*(Angle2**4))/10000 + (273*np.pi**2*(Angle2**5))/500 + 22081/2000
	r2FCR = 14
	r2ECRB = lambda Angle2: (8199*np.pi)/5000 + (44637*Angle2)/2500 - (5073*np.pi*Angle2)/10000 \
							- (471*np.pi*(Angle2**2))/5000 - (28827*(Angle2**2))/10000 - 1407/125
	r2ECRL = lambda Angle2: (74361*np.pi)/10000 + (72089699777459*Angle2)/4398046511104 \
							- (8783*np.pi*Angle2)/5000 + (371*np.pi**2*Angle2)/5000 - (1667*np.pi**2)/1250 \
							- 38517/5000
	r2FCU = 19
	r2FDS = 20
	r2PL = 25
	r2ECU = -23
	r2EDM = -10
	r2EDC = -20
	r2APL = 0

	#R_tranpose Column 3
	r3DELTa = 0
	r3CB = 0
	r3DELTp = 0
	r3BIC = 0
	r3TRI = 0
	r3BRA = 0
	r3BRD = 0
	r3PRO = 0
	r3FCR = lambda Angle3: (3199*Angle3)/2000 + 3301/250
	r3ECRB = lambda Angle3: (21411*Angle3)/10000 - 7562500789275879/562949953421312
	r3ECRL = lambda Angle3: (457*Angle3)/200 - 58583/5000
	r3FCU = lambda Angle3: (13307*(Angle3**2))/10000 + (1869*Angle3)/400 + 1578328710658497/140737488355328
	r3FDS = lambda Angle3: (2099*(Angle3**2))/2000 + (10641*Angle3)/10000 + 5824674283064289/562949953421312
	r3PL = lambda Angle3: (5011*(Angle3**2))/10000 + (13821*Angle3)/10000 + 3749/400
	r3ECU = lambda Angle3: (3883*Angle3)/1250 - 21289/2500
	r3EDM = lambda Angle3: (7603*Angle3)/2500 - 7791/625
	r3EDC = lambda Angle3: (693*Angle3)/400 - 35319/2500
	r3APL = lambda Angle3: 1171/2000 - (171*(Angle3**2))/2000 - (73*Angle3)/2500

	OptimalMuscleLength = np.array([(9.8),     (9.3),      (13.7),      (11.6),      (13.4),    (8.6),   \
                      				(17.3),    (4.9),       (6.3),       (5.9),       (8.1),    (5.1),   \
                       				(8.4),     (6.4),       (6.2),       (6.8),        (7.),    (7.1)] )
	def muscle_1_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1DELTa  									\
					- Angle2Spline.pp_deriv(Time)*r2DELTa  								\
					-Angle3Spline.pp_deriv(Time)*r3DELTa)  								\
         	        /(10*OptimalMuscleLength[0])
		return(velocity)
	def muscle_2_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1CB  								\
                    -Angle2Spline.pp_deriv(Time)*r2CB  									\
                    -Angle3Spline.pp_deriv(Time)*r3CB)  								\
                    /(10*OptimalMuscleLength[1])
		return(velocity)
	def muscle_3_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1DELTp  								\
                    -Angle2Spline.pp_deriv(Time)*r2DELTp  								\
                    -Angle3Spline.pp_deriv(Time)*r3DELTp)  								\
                    /(10*OptimalMuscleLength[2])
		return(velocity)
	def muscle_4_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1BIC  								\
                    -Angle2Spline.pp_deriv(Time)*r2BIC(Angle2Spline.pp_func(Time))  	\
                    -Angle3Spline.pp_deriv(Time)*r3BIC)  								\
                    /(10*OptimalMuscleLength[3])
		return(velocity)
	def muscle_5_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1TRI  								\
                    -Angle2Spline.pp_deriv(Time)*r2TRI(Angle2Spline.pp_func(Time))  	\
                    -Angle3Spline.pp_deriv(Time)*r3TRI)  								\
                    /(10*OptimalMuscleLength[4])
		return(velocity)
	def muscle_6_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1BRA  								\
                    -Angle2Spline.pp_deriv(Time)*r2BRA(Angle2Spline.pp_func(Time))  	\
                    -Angle3Spline.pp_deriv(Time)*r3BRA)  								\
                    /(10*OptimalMuscleLength[5])
		return(velocity)
	def muscle_7_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1BRD  								\
                    -Angle2Spline.pp_deriv(Time)*r2BRD(Angle2Spline.pp_func(Time))  	\
                    -Angle3Spline.pp_deriv(Time)*r3BRD)  								\
                    /(10*OptimalMuscleLength[6])
		return(velocity)
	def muscle_8_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1PRO  								\
                    -Angle2Spline.pp_deriv(Time)*r2PRO(Angle2Spline.pp_func(Time))  	\
                    -Angle3Spline.pp_deriv(Time)*r3PRO)  								\
                    /(10*OptimalMuscleLength[7])
		return(velocity)
	def muscle_9_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1FCR  								\
                    -Angle2Spline.pp_deriv(Time)*r2FCR  								\
                    -Angle3Spline.pp_deriv(Time)*r3FCR(Angle3Spline.pp_func(Time)))  	\
                    /(10*OptimalMuscleLength[8])
		return(velocity)
	def muscle_10_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1ECRB  								\
                    -Angle2Spline.pp_deriv(Time)*r2ECRB(Angle2Spline.pp_func(Time))  	\
                    -Angle3Spline.pp_deriv(Time)*r3ECRB(Angle3Spline.pp_func(Time)))  	\
                    /(10*OptimalMuscleLength[9])
		return(velocity)
	def muscle_11_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1ECRL  								\
                    -Angle2Spline.pp_deriv(Time)*r2ECRL(Angle2Spline.pp_func(Time))  	\
                    -Angle3Spline.pp_deriv(Time)*r3ECRL(Angle3Spline.pp_func(Time))) 	\
                    /(10*OptimalMuscleLength[10])
		return(velocity)
	def muscle_12_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1FCU  								\
                    -Angle2Spline.pp_deriv(Time)*r2FCU  								\
                    -Angle3Spline.pp_deriv(Time)*r3FCU(Angle3Spline.pp_func(Time))) 	\
                    /(10*OptimalMuscleLength[11])
		return(velocity)
	def muscle_13_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1FDS  								\
                    -Angle2Spline.pp_deriv(Time)*r2FDS  								\
                    -Angle3Spline.pp_deriv(Time)*r3FDS(Angle3Spline.pp_func(Time))) 	\
                    /(10*OptimalMuscleLength[12])
		return(velocity)
	def muscle_14_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1PL  									\
                    -Angle2Spline.pp_deriv(Time)*r2PL  									\
                    -Angle3Spline.pp_deriv(Time)*r3PL(Angle3Spline.pp_func(Time))) 		\
                    /(10*OptimalMuscleLength[13])
		return(velocity)
	def muscle_15_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1ECU  								\
                    -Angle2Spline.pp_deriv(Time)*r2ECU  								\
                    -Angle3Spline.pp_deriv(Time)*r3ECU(Angle3Spline.pp_func(Time))) 	\
                    /(10*OptimalMuscleLength[14])
		return(velocity)
	def muscle_16_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1EDM  								\
                    -Angle2Spline.pp_deriv(Time)*r2EDM  								\
                    -Angle3Spline.pp_deriv(Time)*r3EDM(Angle3Spline.pp_func(Time))) 	\
                    /(10*OptimalMuscleLength[15])
		return(velocity)
	def muscle_17_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1EDC  								\
                    -Angle2Spline.pp_deriv(Time)*r2EDC  								\
                    -Angle3Spline.pp_deriv(Time)*r3EDC(Angle3Spline.pp_func(Time))) 	\
                    /(10*OptimalMuscleLength[16])
		return(velocity)
	def muscle_18_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time):
		velocity = (-Angle1Spline.pp_deriv(Time)*r1APL  								\
                    -Angle2Spline.pp_deriv(Time)*r2APL  								\
                    -Angle3Spline.pp_deriv(Time)*r3APL(Angle3Spline.pp_func(Time))) 	\
                    /(10*OptimalMuscleLength[17])
		return(velocity)	
	NormalizedMuscleVelocity = np.array([ 	muscle_1_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_2_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_3_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_4_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_5_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_6_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_7_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_8_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_9_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_10_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_11_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_12_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_13_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_14_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_15_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_16_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_17_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time), \
											muscle_18_velocity(Angle1Spline,Angle2Spline,Angle3Spline,Time) 	],	\
			                                float)
	return(NormalizedMuscleVelocity)
def plot_normalized_muscle_velocity(Angle1Splines,Angle2Splines,Angle3Splines,Time):
	NormalizedMuscleVelocity = normalized_muscle_velocity(Angle1Splines,Angle2Splines,Angle3Splines,Time)
	plt.figure()
	for i in range(18):
		plt.plot(Time,NormalizedMuscleVelocity[i])
	plt.show()
def max_contraction(Angle1Splines,Angle2Splines,Angle3Splines,Time):
	NormalizedMuscleVelocity = normalized_muscle_velocity(Angle1Splines,Angle2Splines,Angle3Splines,Time)
	Zeros = [[0]*len(NormalizedMuscleVelocity[0])]*len(NormalizedMuscleVelocity)
	EccentricContractions = (NormalizedMuscleVelocity>Zeros)*NormalizedMuscleVelocity
	ConcentricContractions = (NormalizedMuscleVelocity<Zeros)*NormalizedMuscleVelocity
	MaxEccContract = np.max(EccentricContractions[:,:5200],1)
	MaxConcContract = np.min(ConcentricContractions[:,:5200],1)
	return(MaxEccContract,MaxConcContract)
def sum_of_squares(Array):
	result = np.sum(Array**2)
	return(result)
def total_sum_of_squares(Angle1Splines,Angle2Splines,Angle3Splines,Time):
	Ecc_Sum_of_Squares = []
	Conc_Sum_of_Squares = []
	for i in range(Angle1Splines.size):
		MaxEccContract,MaxConcContract = max_contraction(Angle1Splines[i],Angle2Splines[i],Angle3Splines[i],Time)
		Ecc_Sum_of_Squares.append(sum_of_squares(MaxEccContract))
		Conc_Sum_of_Squares.append(sum_of_squares(MaxConcContract))
	return(Ecc_Sum_of_Squares,Conc_Sum_of_Squares)
Ecc_Sum_of_Squares,Conc_Sum_of_Squares = total_sum_of_squares(Angle1Splines,Angle2Splines,Angle3Splines,Time)
plt.figure()
plt.hist2d(Ecc_Sum_of_Squares, Conc_Sum_of_Squares, bins=25)
plt.colorbar()
plt.show()
#plot_spline_results(Time,Angle1Spline)
#plot_spline_results(Time,Angle2Spline)
#plot_spline_results(Time,Angle3Spline)


