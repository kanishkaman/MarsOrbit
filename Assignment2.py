import datetime
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Loading and storing data
marsData = pd.read_csv('mars_data.csv')  # Replace with your OWN filepath 
marsData.head(10)

marsHeliocentric_longitude = marsData.values[:, 5:9]
print(marsHeliocentric_longitude)

marsHeliocentric_longitude_InDegree = np.array(marsData['ZodiacIndex'] * 30 + \
                                          marsData['Degree'] + \
                                          marsData['Minute.1']/60.0 + \
                                          marsData['Second']/3600.0)

marsHeliocentric_longitude_InDegree_InRad = marsHeliocentric_longitude_InDegree * math.pi / 180.0

# TIMES
times = list([0])       # We will create a list which will hold all the times
for i in range(1, len(marsData)):
        date1 = datetime.datetime(
            marsData['Year'][i-1],marsData['Month'][i-1],
            marsData['Day'][i-1],marsData['Hour'][i-1],
            marsData['Minute'][i-1])
    
        date2 = datetime.datetime(
            marsData['Year'][i],marsData['Month'][i],
            marsData['Day'][i],marsData['Hour'][i],
            marsData['Minute'][i])
        
        duration = date2 - date1
        numOfDays = duration.days + duration.seconds / (60*60*24)
        times.append(numOfDays)
times = np.array(times)
print(times)

marsHeliocentric_longitude_InDegree = np.array(marsHeliocentric_longitude_InDegree)
oppositions = np.stack((times,marsHeliocentric_longitude_InDegree), axis = 1)
print('Oppositions: \n', oppositions)


# Plot of spokes w.r.t. sun-aries axis 
plt.figure(figsize=(6, 6), dpi=120)

# Loop through the oppositions to plot spokes
for j in range(0, 12):
    x_position = np.cos(np.radians(oppositions[j][1]))
    y_position = np.sin(np.radians(oppositions[j][1]))

    line_x = [0, x_position]
    line_y = [0, y_position]

    plt.plot(line_x, line_y, linewidth=2, linestyle='--')

plt.axis('equal')
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.title('Spokes with Respect to Sun-Aries Axis', fontsize=14)
plt.xlabel('X-axis', fontsize=12)
plt.ylabel('Y-axis', fontsize=12)
plt.show()

# PLOTS TO SUPPORT VISUALISATION
c = 60  # Angle of center from sun-Aries axis (assumed)
orbital_period = 687  # Orbital period of Mars in days
angular_speed = 360 / orbital_period

plt.figure(figsize=(6, 6), dpi=100)

total_time = 0

# Loop through the oppositions to plot spokes
for index in range(1, 12):
    # Update total time based on angular speed
    total_time = (oppositions[index][0] * angular_speed + total_time) % 360
    
    x_position = math.cos(math.radians(total_time))
    y_position = math.sin(math.radians(total_time))

    line_x = [0, x_position]
    line_y = [0, y_position]

    plt.plot(line_x, line_y, linestyle='--', linewidth=2)

plt.axis('equal')
plt.grid(color='gray', linestyle='--', linewidth=0.5)

# Set title and labels
plt.title('Mars Position Over Time', fontsize=14)
plt.xlabel('X-axis', fontsize=12)
plt.ylabel('Y-axis', fontsize=12)
plt.show()

def getIntersectionPoint(h,k,theta, r,c):
        """
        Returns intersection point of 0-centered circle and line from equant\n",
        
        h     : x-coordinate of point through which line passes through
        k     : y-coordinate of point through which line passes through
        theta : angle the line makes with x-axis
        r     : radius of the circle centered at the (1,c) 
        c     : the angle centre makes with sun-aries line
        ell_value : Running variable, different values pf l will give
                    different points on line which intersects with orbit.
        """
        cosineValueOfTheta = math.cos(math.radians(theta))
        sineValueOfTheta = math.sin(math.radians(theta))

        cosineValueOfC = math.cos(math.radians(c))
        sineValueOfC = math.sin(math.radians(c))

        b = 2 * ((h * cosineValueOfTheta) + (k * sineValueOfTheta)
                -(cosineValueOfC * cosineValueOfTheta)
                -(sineValueOfC * sineValueOfTheta))
        
        c1 = h**2 + k**2 + 1 - (2 * h * cosineValueOfC) - (2 * k * sineValueOfC) - r**2
        l1 = -b / 2
        
        try:
            l2 = math.sqrt(b ** 2 -(4 * c1)) / 2
        except :
            l2 = 0
            
        root1 = l1 + l2
        root2 = l1 - l2
        
        if root1 > 0:
            ell_value = root1
        else:
            ell_value = root2
        
        return (h + ell_value * cosineValueOfTheta), (k + ell_value * sineValueOfTheta)

# QUESTION 1: MARS EQUANT MODEL
def MarsEquantModel(c,r,e1,e2,z,s,oppositions):
        """
        Return 12 errors of angle delta wrt to 12 oppositions and max error among these 12
        """
        errors = []
        xpos = list()
        ypos = list()
        
        h = e1 * math.cos(math.radians(e2 + z))
        k = e1 * math.sin(math.radians(e2 + z))
    
        thetaNew = z
        for i in range(12):
            theta = (s * times[i]) + thetaNew
            x,y = getIntersectionPoint(h,k,theta,r,c)
            xpos.append(x)
            ypos.append(y)
            angle = np.degrees(np.arctan2(y,x))%360
            # print("projection of mars -in degree for spoke ",i,":-",angle,"\n")
            errors.append(abs(oppositions[i][1] - angle))
            thetaNew = theta
        maxError = max(errors)
        return errors, maxError

errors, maxError = MarsEquantModel(149,8.599999,1.60000000001,93.2,55.800000001,0.524,oppositions)
print("12 errors-:",errors)
print("----------------------------------------------------------------------------------------------------------------")
print('maximum error:-',maxError)

# QUESTION 2: BEST ORBIT INNER PARAMS
def bestOrbitInnerParams(r,s,oppositions):          
        maxError = 1e15
        for c in (np.arange(149,149.4,0.01)):
            for e2 in np.arange(93,94,0.1):
                for z in np.arange(55.5,56,0.1):
                    for e1 in np.arange(1.45, 1.6,0.01):
                        errors, max_Error = MarsEquantModel(c,r,e1,e2,z,s,oppositions)
                        if(maxError > max_Error):
                            maxError = round(max_Error,4)
                            max_c = c
                            max_e1 = round(e1,4)
                            max_e2 = e2
                            max_z = z
                            ErrorList = errors
                            print("C value:-",max_c," E1 value:-",max_e1," E2 value:-",max_e2," z value:-",max_z," Maximum Error:-",maxError)
        return max_c,max_e1,max_e2,max_z,ErrorList,maxError

c,e1,e2,z,errors,maxError = bestOrbitInnerParams(9,0.524,oppositions)
print("----------------------------------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------------------------------")
print("C value:-",c," E1 value:-",e1," E2 value:-",e2," z value:-",z," Error List:-",errors," Maximum Error:-",maxError)

# QUESTION 3: BEST S
def bestS(r, oppositions):
        leastError = 1e15
        for s in np.arange(680, 689,0.2):
            max_c,max_e1,max_e2,max_z,errors,maxError= bestOrbitInnerParams(r,360/s,oppositions)
            if(leastError > maxError):
                bestS = 360/s
                ErrorList = errors
                leastError = round(maxError,4)
                print("S value:-",bestS,"C value:-",max_c," E1 value:-",max_e1," E2 value:-",max_e2," z value:-",max_z," Maximum Error:-",maxError)
        return bestS,ErrorList,leastError

s,errors,maxError = bestS(8,oppositions)
print("----------------------------------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------------------------------")
print("S value:-",s," Error List:-",errors," Maximum Error:-",maxError)

# QUESTION 4: BEST R
def bestR(s, oppositions):
        leastError = 1e15
        for r in np.arange(5, 10,0.1):
            max_c,max_e1,max_e2,max_z,errors,maxError= bestOrbitInnerParams(r,s,oppositions)
            if(leastError > maxError):
                bestR = r
                ErrorList = errors
                leastError = round(maxError,4)
                print("R value:-",bestR,"C value:-",max_c," E1 value:-",max_e1," E2 value:-",max_e2," z value:-",max_z," Maximum Error:-",maxError),
        return bestR,ErrorList,leastError

r,errors,maxError = bestR(0.524,oppositions)
print("----------------------------------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------------------------------")
print("Best R value:- ",r," Error List:- ",errors," Max Error:- ",maxError)

# QUESTION 5: BEST MARS ORBIT PARAMS
def bestMarsOrbitParams(oppositions):
    bestError = 1e15
    for r in (np.arange(8,8.3,0.01)):
        for s in np.arange(686.5,687,0.1):  # 686.93 days is the exact value.
            max_c,max_e1,max_e2,max_z,errors,maxError= bestOrbitInnerParams(r,360/s,oppositions)
            if(bestError > maxError):
                optimumCVal = max_c
                optimume1Val = max_e1
                optimume2Val = max_e2
                optimumZVal = max_z
                bestErrorList= errors
                bestError = maxError
                optimumRVal = r
                optimumSVal = 360/s
                print("Best Error:- ",bestError,"optimumCVal:- ",optimumCVal,"optimume1Val:- ",optimume1Val,"optimume2Val:- ",optimume2Val,"optimumZVal",optimumZVal,"optimumRVal:- ",optimumRVal,"optimumSVal",optimumSVal)
    return optimumRVal,optimumSVal,optimumCVal,optimume1Val,optimume2Val,optimumZVal,bestErrorList,bestError

r,s,c,e1,e2,z,errors,maxError = bestMarsOrbitParams(oppositions)
print("----------------------------------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------------------------------")
print("Best Error: ",maxError,"optimumCVal: ",c,"optimume1Val: ",e1,"optimume2Val: ",e2,"optimumZVal: ",z,"optimumRVal: ",r,"optimumSVal: ",s)

# DISPLAYING RESULTS IN A TABULAR FORMAT
dataa = {
    "Parameter": ["Best Error", "Optimum C Value", "Optimum e1 Value", "Optimum e2 Value", "Optimum Z Value", "Optimum R Value", "Optimum S Value"],
    "Value": [maxError, c, e1, e2, z, r, s]
}
daf = pd.DataFrame(dataa)
print(daf)


# FUNCTION TO PLOT MARS ORBIT
def plotMarsOrbit(c, r, e1, e2, z, oppositions):
    """
        Plot Mars orbit given the center of the circle and radius along with spokes from equant and sun.
    """

    fig, ax = plt.subplots(figsize=(12, 12))

    ax.set_xlim((-12, 12))
    ax.set_ylim((-12, 12))

    CentreXPos = math.cos(math.radians(c))
    CentreYPos = math.sin(math.radians(c))
    
    # Calculate equant position
    h = e1 * math.cos(math.radians(e2 + z))
    k = e1 * math.sin(math.radians(e2 + z))
    
    xpos = []
    ypos = []
    
    thetaNew = z
    for i in range(12):
        theta = (s * oppositions[i][0]) + thetaNew
        x, y = getIntersectionPoint(h, k, theta, r, c)
        plt.scatter(x, y, color='black')
        plt.plot([h, x], [k, y], linestyle='dotted')
        xpos.append(x)
        ypos.append(y)
        thetaNew = theta
        
    # Create the orbit circle
    orbit = plt.Circle((CentreXPos, CentreYPos), r, color='cyan', fill=False)
    
    for i in range(12):
        xpos = math.cos(math.radians(oppositions[i][1]))
        ypos = math.sin(math.radians(oppositions[i][1]))
        
        x = [0, xpos * 25]  # Increased length of the line
        y = [0, ypos * 25]
 
        plt.plot(x, y)
    
    # Annotate the center, sun, and equant
    ax.annotate('Centre', xy=(CentreXPos, CentreYPos), fontsize=15, color='blue', fontweight='bold')
    ax.annotate('Sun', xy=(0, 0), fontsize=15, color='red', fontweight='bold')
    ax.annotate('Equant', xy=(h, k), fontsize=15, color='green', fontweight='bold')

    ax.scatter([0, CentreXPos, h], [0, CentreYPos, k], color='purple')
    
    # Add the orbit to the plot
    ax.add_artist(orbit)
    plt.title('Predicted Orbit of Mars', fontsize=18)
    plt.grid(True)
    plt.show()

plotMarsOrbit(149, 8.599999, 1.6, 93.2, 55.8, oppositions)
