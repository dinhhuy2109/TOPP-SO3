from openravepy import *
from numpy import *

from toppso3.SO3RRT import *
from toppso3 import Utils
from toppso3 import lie

import time
import TOPP
from TOPP import TOPPpy
from TOPP import TOPPbindings
from TOPP import Trajectory
from pylab import *
import scipy.optimize
from mpl_toolkits.mplot3d import Axes3D


ion()



env = Environment()
# This model was downloaded from http://nasa3d.arc.nasa.gov/models/printable
env.Load("../MESSENGER/messengerWithEnv.xml")
env.SetViewer('qtcoin')

robot = env.GetBodies()[0]


phi = pi

R0 = eye(3)
q0 = quatFromRotationMatrix(R0)
omega0 = zeros(3)

q1 = array([cos(phi/2.),0,0,sin(phi/2.)])
omega1 = zeros(3)

taumax = ones(3)
vmax = ones(3)
inertia = eye(3)
################################## BiRRT planner #################################

vertex_beg = Vertex(Config(q0,omega0), FW)
vertex_end = Vertex(Config(q1,omega1), BW)
biRRTinstance = RRTPlanner(vertex_beg, vertex_end, robot)

allottedtime = 600
biRRTinstance.Run(allottedtime)

Rlist = biRRTinstance.GenFinalRotationMatrixList()
Trajlist = biRRTinstance.GenFinalTrajList()
lietraj = lie.LieTraj(Rlist,Trajlist)

ion()

## Visualize 
# M = eye(4)
# for t in linspace(0, lietraj.duration, 1000): 
#     M[:3,:3] = lietraj.EvalRotation(t)
#     robot.SetTransform(M)
#     isincollision = (env.CheckCollision(robot, CollisionReport()))
#     if (isincollision):
#         print "in collision", " ", t, "/" , lietraj.duration
#     time.sleep(0.01)

################################# TOPP #############################################
discrtimestep= 1e-2
constraintsstring = str(discrtimestep)
constraintsstring += "\n" + ' '.join([str(v) for v in taumax])
for v in inertia:
    constraintsstring += "\n" + ' '.join([str(i) for i in v])
# Note that, when Inertia is an Identity matrix, angular accelerations are the same as torques
print "\033[93mRunning TOPP", "\033[0m"

t_topp_start = time.time()
traj = Trajectory.PiecewisePolynomialTrajectory.FromString(Utils.TrajStringFromTrajList(Trajlist))

abc = TOPPbindings.RunComputeSO3Constraints(str(traj),constraintsstring)
a,b,c = lie.Extractabc(abc)
# a,b,c = lie.ComputeSO3Constraints(traj, taumax, discrtimestep) #This is the implementation of computing SO3Constraints in Python
topp_inst = TOPP.QuadraticConstraints(traj, discrtimestep, vmax, list(a), list(b), list(c))

x = topp_inst.solver

ret = x.RunComputeProfiles(0,0)
if ret == 1:
    x.ReparameterizeTrajectory()
    x.WriteResultTrajectory()

traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
t_topp_end = time.time()

print "\033[1;32mRunning time:",t_topp_end-t_topp_start, "sec.\033[0m"
print "\033[93mDone", "\033[0m"
lietraj1 = lie.SplitTraj2(Rlist, traj1)

#---Visualize----
# M = eye(4)

# for t in linspace(0, lietraj1.duration, 1000): 
#     M[:3,:3] = lietraj1.EvalRotation(t)
#     robot.SetTransform(M)
#     isincollision = (env.CheckCollision(robot, CollisionReport()))
#     if (isincollision):
#         print "in collision", " ", t, "/" , lietraj1.duration
#     time.sleep(0.01)

################################ SHORTCUTTING ############################

print "\033[93mRunning SHORTCUTTING", "\033[0m"

taumax = ones(3)
vmax = ones(3)
lietraj2 = Utils.Shortcut(robot, taumax, vmax, lietraj1, 200, -1, 0, -1, inertia)

print "\033[93mDone", "\033[0m"

print "\033[1;94mFinal trajectory duration: ", lietraj2.duration, " sec.\033[0m"

#---Visualize----
M = eye(4)
for t in linspace(0, lietraj2.duration, 1000): 
    M[:3,:3] = lietraj2.EvalRotation(t)
    robot.SetTransform(M)
    isincollision = (env.CheckCollision(robot, CollisionReport()))
    if (isincollision):
        print "in collision", " ", t, "/" , lietraj2.duration
    time.sleep(0.01)

# lietraj2.Plot(0.01,0,vmax,taumax,taumax,inertia)

################# SAVE LIETRAJ #########################################
#Utils.SaveLietrajAsTextFiles(lietraj1, "Rlist1.txt", "trajlist1.txt")
#Utils.SaveLietrajAsTextFiles(lietraj2, "Rlist2.txt", "trajlist2.txt")


#-------------------- Plotting the MVC and the profiles --------------------#
# x.WriteProfilesList()
# x.WriteSwitchPointsList()
# profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
# switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
# TOPPpy.PlotProfiles(profileslist, switchpointslist, 4)


##########################LOAD LIETRAJ #################################
# lietraj4 = Utils.ReadLieTrajFiles("Rlist2.txt", "trajlist2.txt")
# for t in linspace(0, lietraj4.duration, 1000): 
#     M[:3,:3] = lietraj4.EvalRotation(t)
#     robot.SetTransform(M)
#     isincollision = (env.CheckCollision(robot, CollisionReport()))
#     if (isincollision):
#         print "in collision", " ", t, "/" , lietraj4.duration
#     time.sleep(0.01)

