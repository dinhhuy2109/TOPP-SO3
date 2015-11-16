from openravepy import *
from numpy import *

from SE3RRT import *

import Utils

import time
import lie
import TOPP
from TOPP import TOPPpy
from TOPP import Trajectory
from pylab import *
import scipy.optimize
from mpl_toolkits.mplot3d import Axes3D

import pdb

ion()



env = Environment()
# This model was downloaded from http://nasa3d.arc.nasa.gov/models/printable
env.Load("MESSENGER/messengerWithEnvSE3.xml")
env.SetViewer('qtcoin')

robot = env.GetBodies()[0]


phi = pi

R0 = eye(3)
q0 = quatFromRotationMatrix(R0)
omega0 = zeros(3)

q1 = array([cos(phi/2.),0,0,sin(phi/2.)])
omega1 = zeros(3)

taumax = ones(3)
vmax = ones(6)
fmax = ones(3)

tran0 = array([0,0,0])
vtran0 = zeros(3)

tran1 = array([0,0.03,0])
vtran1 = zeros(3)


# pdb.set_trace()
# ################################## BiRRT planner #################################

vertex_beg = Vertex(Config(q0,tran0, omega0, vtran0), FW)
vertex_end = Vertex(Config(q1,tran1, omega1, vtran1), BW)
biRRTinstance = RRTPlanner(vertex_beg, vertex_end, robot)

allottedtime = 600
biRRTinstance.Run(allottedtime)

Rlist = biRRTinstance.GenFinalRotationMatrixList()
TrajRotlist = biRRTinstance.GenFinalTrajList()
lietraj = lie.LieTraj(Rlist,TrajRotlist)

transtraj = Trajectory.PiecewisePolynomialTrajectory.FromString(biRRTinstance.GenFinalTrajTranString())


ion()

#---Visualize----
# M = eye(4)
# for t in linspace(0, lietraj.duration, lietraj.duration*100): 
    # M[:3,:3] = lietraj.EvalRotation(t)
    # M[:3, 3] = transtraj.Eval(t)
    # robot.SetTransform(M)
    # isincollision = (env.CheckCollision(robot, CollisionReport()))
    # if (isincollision):
    #     print "in collision", " ", t, "/" , lietraj.duration
    # time.sleep(0.01)

################################# TOPP ############################################# 
print "\033[93mRunning TOPP", "\033[0m"

t_topp_start = time.time()
rottraj = Trajectory.PiecewisePolynomialTrajectory.FromString(Utils.TrajStringFromTrajList(TrajRotlist))
se3traj = Utils.SE3TrajFromTransandSO3(transtraj, rottraj)

discrtimestep = 1e-2

a,b,c = Utils.ComputeSE3Constraints(se3traj, taumax, fmax, discrtimestep)
topp_inst = TOPP.QuadraticConstraints(se3traj, discrtimestep, vmax, list(a), list(b), list(c))
x = topp_inst.solver
ret = x.RunComputeProfiles(0,0)
if ret == 1:
    x.ReparameterizeTrajectory()
    x.WriteResultTrajectory()

se3traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)


t_topp_end = time.time()

print "\033[1;32mRunning time:",t_topp_end-t_topp_start, "sec.\033[0m"
print "\033[93mDone", "\033[0m"

transtraj1, rottraj1 = Utils.TransRotTrajFromSE3Traj(se3traj1)
lietraj1 = lie.SplitTraj2(Rlist, rottraj1)

#---Visualize----
# M = eye(4)
# for t in linspace(0, lietraj1.duration,5*100): 
#     M[:3,:3] = lietraj1.EvalRotation(t)
#     M[:3, 3] = transtraj1.Eval(t)
#     robot.SetTransform(M)
#     isincollision = (env.CheckCollision(robot, CollisionReport()))
#     if (isincollision):
#         print "in collision", " ", t, "/" , lietraj1.duration
#     time.sleep(0.01)

############################### SHORTCUTTING ################################
print "\033[93mRunning SHORTCUTING", "\033[0m"
se3traj2, Rlist2 = Utils.SE3Shortcut(robot, taumax, fmax, vmax, se3traj1, Rlist, 200)#, -1,0,-1,1)

#print se3traj2.duration

transtraj2, rottraj2 = Utils.TransRotTrajFromSE3Traj(se3traj2)
lietraj2 = lie.SplitTraj2(Rlist2, rottraj2)

M = eye(4)
for t in linspace(0, lietraj2.duration, 5*100): 
    M[:3,:3] = lietraj2.EvalRotation(t)
    M[:3, 3] = transtraj2.Eval(t)
    robot.SetTransform(M)
    isincollision = (env.CheckCollision(robot, CollisionReport()))
    if (isincollision):
        print "in collision", " ", t, "/" , lietraj2.duration
    time.sleep(0.01)


print "\033[1;94mFinal trajectory duration: ", se3traj2.duration, " sec.\033[0m"

# Utils.PlotSE3(se3traj2, Rlist2, 0.01, 0,vmax,taumax,taumax,fmax,np.eye(3))

# #################### SAVE se3traj ############################################
# Utils.SaveSE3trajAsTextFiles(se3traj1, Rlist,  "se3Rlist1.txt", "se3trajlist1.txt")
# Utils.SaveSE3trajAsTextFiles(se3traj2, Rlist2, "se3Rlist2.txt", "se3trajlist2.txt")


###########################LOAD SE3TRAJ ##########################################
# se3traj3, Rlist3 = Utils.ReadSE3TrajFiles("se3Rlist2.txt", "se3trajlist2.txt")

# transtraj3, rottraj3 = Utils.TransRotTrajFromSE3Traj(se3traj3)
# lietraj3 = lie.SplitTraj2(Rlist3, rottraj3)

# M = eye(4)
# for t in linspace(0, lietraj3.duration, 5*100): 
#     M[:3,:3] = lietraj3.EvalRotation(t)
#     M[:3, 3] = transtraj3.Eval(t)
#     robot.SetTransform(M)
#     isincollision = (env.CheckCollision(robot, CollisionReport()))
#     if (isincollision):
#         print "in collision", " ", t, "/" , lietraj3.duration
#     time.sleep(0.01)
