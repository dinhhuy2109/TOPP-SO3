from openravepy import *
from pylab import *

from numpy import *

import lie
import time 

import string
import StringIO

import TOPP
from TOPP import TOPPpy
from TOPP import TOPPbindings
from TOPP import Trajectory
from TOPP import Utilities

import pdb ########################
import matplotlib.pyplot as plt


def QuatDistance(quat0, quat1): 
    rotationweight = 1
    innerProduct = dot(quat0,quat1)
    quatDistance = rotationweight*(1-abs(innerProduct))
    return quatDistance

def SO3Distance(R0, R1): # bi-invariance
    return linalg.norm(lie.logvect(dot(R0.T,R1)))

def R3Distance(b0, b1):
    return linalg.norm(b1-b0)

def SE3Distance(X0, X1, c = None, d = None): # left invariance
    R0 = X0[:3,:3]
    R1 = X1[:3,:3]
    b0 = X0[:3,3]
    b1 = X1[:3,3]
    if (c == None):
        c = 1
    else: c = c
    if (d == None):
        d = 1
    else: d = d
    return sqrt(c*(SO3Distance(R0,R1)**2) + d*(R3Distance(b0,b1)**2))
    
################## interpolate translation ####################################
def TrajString3rdDegree(q_beg, q_end, qs_beg, qs_end, duration):
    trajectorystring = ''
    ndof = len(q_beg)
    
    trajectorystring += "%f\n%d"%(duration, ndof)

    for k in range(ndof):
        a, b, c, d = Utilities.Interpolate3rdDegree(q_beg[k], q_end[k], qs_beg[k], qs_end[k], duration)
        trajectorystring += "\n%f %f %f %f"%(d, c, b, a)
    return trajectorystring

#################### SE3 traj ##################################################

def SE3TrajFromTransandSO3(transtraj, rtraj): # same chunk.duration
    #return duration-dimension-trans polynomial- rot polynomial
    if len(transtraj.chunkslist) != len(rtraj.chunkslist):
        print 'error'
        return 0
    clist = []
    for c in transtraj.chunkslist:
        plist = []
        for i in range(3):
            plist.append(c.polynomialsvector[i])
        for i in range(3):
            rc = rtraj.chunkslist[len(clist)]
            plist.append(rc.polynomialsvector[i])
        chunk = Trajectory.Chunk(c.duration, plist)
        clist.append(chunk)
    return Trajectory.PiecewisePolynomialTrajectory(clist)

###################### Decompose SE3 traj to ROT and Trans traj ################
def TransRotTrajFromSE3Traj(SE3traj):
    transclist = []
    rclist = []
    for c in SE3traj.chunkslist:
        transchunk = Trajectory.Chunk(c.duration, c.polynomialsvector[:3])
        transclist.append(transchunk)
        rchunk = Trajectory.Chunk(c.duration, c.polynomialsvector[3:])
        rclist.append(rchunk)
    transtraj = Trajectory.PiecewisePolynomialTrajectory(transclist)
    rtraj = Trajectory.PiecewisePolynomialTrajectory(rclist)
    return transtraj, rtraj

##########################SE3 constraint ########################################
def ComputeSE3Constraints(SE3traj, taumax, fmax, discrtimestep, I = None, m = None):
    ndiscrsteps = int((SE3traj.duration + 1e-10) / discrtimestep) + 1
    a = zeros((ndiscrsteps,12))
    b = zeros((ndiscrsteps,12))
    c = zeros((ndiscrsteps,12))
    transtraj, rtraj = TransRotTrajFromSE3Traj(SE3traj)
    for i in range(ndiscrsteps):
        #rotconstraints
        t = i * discrtimestep
        r = rtraj.Eval(t)
        rd = rtraj.Evald(t)
        rdd = rtraj.Evaldd(t)
        nr = linalg.norm(r)
        nr2 = nr*nr
        nr3 = nr2*nr
        nr4 = nr3*nr
        nr5 = nr4*nr
        R = lie.skewfromvect(r)

        snr = sin(nr)
        cnr = cos(nr)
        rcrd = cross(r,rd)
        rdrd = dot(r,rd)

        Amat =  eye(3) - (1-cnr)/nr2*R + (nr-snr)/nr3*dot(R,R)
        C1 = (nr-snr)/nr3 * cross(rd,rcrd)
        C2 = -(2*cnr+nr*snr-2)/nr4 * rdrd*rcrd
        C3 = (3*snr-nr*cnr - 2*nr)/nr5 * rdrd*cross(r,rcrd)
        C = C1+C2+C3

        Ard = dot(Amat,rd)
        if I is None:            
            at = Ard
            bt = dot(Amat,rdd) + C
        else:
            at = dot(I,Ard)
            bt = dot(I,dot(Amat,rdd)) + dot(I,C) + cross(Ard,dot(I,Ard))
        
        a[i,3:6] = at
        a[i,9:12] = -at
        b[i,3:6] = bt
        b[i,9:12] = -bt
        c[i,3:6] = -taumax
        c[i,9:12] = -taumax

        #transconstraints
        td = transtraj.Evald(t)
        tdd = transtraj.Evaldd(t)
        if m is None:
            at = td
            bt = tdd

        a[i,:3] = at
        a[i,6:9] = -at
        b[i,:3] = bt
        b[i,6:9] = -bt
        c[i,:3] = -fmax
        c[i,6:9] = -fmax
    return a, b, c

######################## se3 traj collision checking ########################

def CheckCollisionSE3Traj( robot, transtraj, rtraj, R_beg,  checkcollisiontimestep = 1e-3):
    """CheckCollisionSE3Traj accepts a robot and trans, rot trajectory object as its inputs.
       (checkcollisiontimestep is set to 1e-3 as a default value)
       It returns True if any config along the traj is IN-COLLISION.
    """
    env = robot.GetEnv()
    for s in np.arange(0, transtraj.duration, checkcollisiontimestep):
        with robot:
            transformation = eye(4)
            transformation[0:3,0:3] = lie.EvalRotation(R_beg, rtraj, s)
            transformation[0:3,3] = transtraj.Eval(s)
            robot.SetTransform(transformation)           
            isincollision = (env.CheckCollision(robot, CollisionReport()))
            #print  "s =", s, " ", isincollision
            if (isincollision):
                return True
    with robot:
        robot.SetTransform(transformation)
        isincollision = (env.CheckCollision(robot, CollisionReport()))
        if (isincollision):
            return True
        else:
            return False


######################### SE3 shortcutting ##################################
def SE3Shortcut(robot, taumax, fmax, vmax, se3traj, Rlist, maxiter, expectedduration = -1,  meanduration = 0, upperlimit = -1, plotdura = None):
    if plotdura == 1:
        plt.axis([0, maxiter, 0, se3traj.duration])
        plt.ion()
        plt.show()
        ylabel('Trajectory duration (s)')
        xlabel('Iteration')

    t_sc_start = time.time()
    originalduration =  se3traj.duration
    #return shortcuted traj
    if upperlimit < 0:
        dur = se3traj.duration
        upperlimit = se3traj.duration
    else:
        dur = upperlimit
    attempt = 0

    ## for shortcutting
    integrationtimestep = 1e-2             
    reparamtimestep = 1e-2                  
    passswitchpointnsteps = 5                
    discrtimestep = 1e-2                    
    assert(dur > 10.0*discrtimestep)
    
    ncollision = 0
    nnotretimable = 0 
    nnotshorter = 0
    
    transtraj, rtraj = TransRotTrajFromSE3Traj(se3traj)
    lietraj = lie.SplitTraj(Rlist, rtraj)
   

    for it in range(maxiter):
        if plotdura == 1:
            plt.scatter(it, se3traj.duration)
            plt.draw()
        #transtraj, rtraj = TransRotTrajFromSE3Traj(se3traj)
        #lietraj = lie.SplitTraj2(Rlist, rtraj)
        if (expectedduration > 0): # check, if newlietraj.duration is short enough, stop SHORTCUTING
            if (se3traj.duration < expectedduration):
                print "\033[1;32mTrajectory's duration is already shorter than expected time --> stop shortcuting\033[0m"
                break
        if (dur < discrtimestep):
            print "[Utils::Shortcut] trajectory duration is less than discrtimestep.\n"
            break ## otherwise, this will cause an error in TOPP        
        
        ## select an interval for shortcutting
        t0 = random.rand()* dur
        
        if meanduration == 0:
            meanduration = dur - t0
            
        T = random.rand()*min(meanduration,dur - t0)
        t1 = t0 + T

        while (T < 2.0*discrtimestep):
            t0 = random.rand()*dur
            if meanduration == 0:
                meanduration = dur - t0
                
            T = random.rand()*min(meanduration, dur - t0)
            t1 = t0 + T

            if t1 > upperlimit:
                t1 = upperlimit
                if (t1 < t0):
                    temp = t0
                    t0 = t1
                    t1 = temp
                    T = t1 - t0

        # print "\n\nShortcutting iteration", it + 1
        # print t0, t1, t1- t0       
        # interpolate from t0 to t1
        R_beg = lietraj.EvalRotation(t0)
        R_end = lietraj.EvalRotation(t1)
        omega0 = lietraj.EvalOmega(t0)
        omega1 = lietraj.EvalOmega(t1)
        shortcutrtraj = lie.InterpolateSO3(R_beg,R_end,omega0,omega1, T)

        t_beg = transtraj.Eval(t0)
        t_end = transtraj.Eval(t1)
        v_beg = transtraj.Evald(t0)
        v_end = transtraj.Evald(t1)
        
        shortcuttranstraj = Trajectory.PiecewisePolynomialTrajectory.FromString(TrajString3rdDegree(t_beg,t_end,v_beg,v_end, T))
        
        shortcutse3traj = SE3TrajFromTransandSO3(shortcuttranstraj, shortcutrtraj)
        #check feasibility only for the new portion
        
        isincollision = CheckCollisionSE3Traj(robot, shortcuttranstraj, shortcutrtraj, R_beg, discrtimestep)
        if (not isincollision):
            a,b,c = ComputeSE3Constraints(shortcutse3traj, taumax, fmax, discrtimestep)
            topp_inst = TOPP.QuadraticConstraints(shortcutse3traj, discrtimestep, vmax, list(a), list(b), list(c))
            x = topp_inst.solver
            ret = x.RunComputeProfiles(1,1) 
            if (ret == 1):
                x.resduration
                ## check whether the new one has shorter duration
                if (x.resduration + 0.1 < T): #skip if not shorter than 0.1 s
                    
                    x.ReparameterizeTrajectory()
                    x.WriteResultTrajectory()
                    TOPPed_shortcutse3traj = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                    TOPPed_shortcuttranstraj, TOPPed_shortcutrtraj = TransRotTrajFromSE3Traj(TOPPed_shortcutse3traj)
                    newlietraj = ReplaceTrajectorySegment(lietraj, TOPPed_shortcutrtraj , t0, t1)
                   
                    newtranstraj = ReplaceTransTrajectorySegment(transtraj, TOPPed_shortcuttranstraj, t0,t1)

                    #####################################################
                    newrtraj = Trajectory.PiecewisePolynomialTrajectory.FromString(TrajStringFromTrajList(newlietraj.trajlist))
                    newse3traj = SE3TrajFromTransandSO3(newtranstraj,newrtraj)

                    Rlist = newlietraj.Rlist
                    rtraj = newrtraj
                    transtraj = newtranstraj
                    lietraj = newlietraj
                    se3traj = newse3traj
                    
                    dur = se3traj.duration

                    #print "*******************************************"
                    print "Success at iteration",it + 1,":", t0, t1,"Deta_t:", t1 - t0 - x.resduration
                    attempt += 1
                    #print "T:", nnotretimable, "; S:", nnotshorter , "; C:", ncollision , "; OK:", attempt
                    #print "*******************************************"
                else:
                    # print "Not shorter"
                    nnotshorter += 1
            else: 
                # print "Not retimable"
                nnotretimable += 1
        else:
            # print "Collision"
            ncollision += 1

            # print "T:", nnotretimable, "; S:", nnotshorter , "; C:", ncollision , "; OK:", attempt
    print "\033[1;32mT:", nnotretimable, "; S:", nnotshorter , "; C:", ncollision , "; OK:", attempt, "\033[0m"
    print "\033[1;32m", originalduration - se3traj.duration ,"sec. shorter\033[0m"
    t_sc_end = time.time()
    print "\033[1;32mRunning time:",t_sc_end-t_sc_start, "sec.\033[0m"
    
    return se3traj, Rlist



#############################

def ReplaceTransTrajectorySegment(originaltranstraj, transtrajsegment, t0,t1):
    """ReplaceTransTrajectorySegment replaces the segment (t0, t1) in the (arbitrary degree) originaltranstraj 
    with an (arbitrary degree) transtrajsegment.
    """
    assert(t1 > t0)
    
    newchunkslist = []
    i0, rem0 = originaltranstraj.FindChunkIndex(t0)
    i1, rem1 = originaltranstraj.FindChunkIndex(t1)
             
    ## check if t0 falls in the first chunk. 
    ## if not, insert chunk 0 to chunk i0 - 1 into newchunkslist
    if i0 > 0:
        for c in originaltranstraj.chunkslist[0: i0]:
            newchunkslist.append(c)

    ## remainderchunk0
    remchunk0 = Trajectory.Chunk(rem0, originaltranstraj.chunkslist[i0].polynomialsvector)
    newchunkslist.append(remchunk0)

    ## insert transtrajsegment
    for c in transtrajsegment.chunkslist:
        newchunkslist.append(c)

    ## remainderchunk1
    newpoly_list = []
    for p in originaltranstraj.chunkslist[i1].polynomialsvector:
        ## perform variable changing of p(x) = a_n(x)^n + a_(n-1)(x)^(n-1) + ...
        ## by x = y + rem1
        
        a = p.q ## coefficient vector with python convention (highest degree first)
        ## a is a poly1d object
        r = a.r ## polynomial roots
        for i in range(len(r)):
            r[i] = r[i] - rem1
        b = np.poly1d(r, True) ## reconstruct a new polynomial from roots
        ## b is a poly1d object
        b = b*a.coeffs[0] ## multiply back by a_n *** this multiplication does not commute
        
        newpoly = Trajectory.Polynomial(b.coeffs.tolist()[::-1]) ## TOPP convention is weak-term-first
        newpoly_list.append(newpoly)
    remchunk1 = Trajectory.Chunk(originaltranstraj.chunkslist[i1].duration - rem1, newpoly_list)
    newchunkslist.append(remchunk1)
    
    ## insert remaining chunks
    if i1 < len(originaltranstraj.chunkslist) - 1:
        for c in originaltranstraj.chunkslist[i1 + 1: len(originaltranstraj.chunkslist)]:
            newchunkslist.append(c)

    return Trajectory.PiecewisePolynomialTrajectory(newchunkslist)




############################# traj collision checking ###############################

def CheckCollisionTraj(robot, trajectory, R_beg, checkcollisiontimestep = 1e-3):
    """CheckCollisionTraj accepts a robot and a trajectory object as its inputs.
       (checkcollisiontimestep is set to 1e-3 as a default value)
       It returns True if any config along the traj is IN-COLLISION.
    """
    env = robot.GetEnv()
    traj = trajectory
    for s in np.arange(0, traj.duration, checkcollisiontimestep):
        with robot:
            transformation = eye(4)
            transformation[0:3,0:3] = lie.EvalRotation(R_beg, traj, s)
            robot.SetTransform(transformation)           
            isincollision = (env.CheckCollision(robot, CollisionReport()))
            #print  "s =", s, " ", isincollision
            if (isincollision):
                return True
    with robot:
        robot.SetTransform(transformation)
        isincollision = (env.CheckCollision(robot, CollisionReport()))
        if (isincollision):
            return True
        else:
            return False



############################# SHORTCUTING SO3 ############################
def Shortcut(robot, taumax, vmax, lietraj,  maxiter, expectedduration = -1,  meanduration = 0, upperlimit = -1, inertia = None, trackingplot = None):
    if trackingplot == 1:
        plt.axis([0, maxiter, 0, lietraj.duration])
        plt.ion()
        plt.show()
        ylabel('Trajectory duration (s)')
        xlabel('Iteration')
    
    
    t_sc_start = time.time()
    originalduration =  lietraj.duration
    #return shortcuted traj
    if upperlimit < 0:
        dur = lietraj.duration
        upperlimit = lietraj.duration
    else:
        dur = upperlimit
    attempt = 0

    ## for shortcutting
    integrationtimestep = 1e-2            
    reparamtimestep = 1e-2                
    passswitchpointnsteps = 5            
    discrtimestep = 1e-2                 

    constraintsstring = str(discrtimestep)
    constraintsstring += "\n" + ' '.join([str(v) for v in taumax])
    if not(inertia is None):
        for v in inertia:
            constraintsstring += "\n" + ' '.join([str(i) for i in v])

    assert(dur > 10.0*discrtimestep)
    

    ncollision = 0
    nnotretimable = 0 
    nnotshorter = 0

    for it in range(maxiter):
        if trackingplot == 1:
            plt.scatter(it, lietraj.duration)
            plt.draw()
        if (expectedduration > 0): # check, if newlietraj.duration is short enough, stop SHORTCUTING
            if (lietraj.duration < expectedduration):
                print "\033[1;32mTrajectory's duration is already shorter than expected time --> stop shortcuting\033[0m"
                break
        if (dur < discrtimestep):
            print "[Utils::Shortcut] trajectory duration is less than discrtimestep.\n"
            break ## otherwise, this will cause an error in TOPP        
        
        ## select an interval for shortcutting
        t0 = random.rand()* dur
        
        if meanduration == 0:
            meanduration = dur - t0
            
        T = random.rand()*min(meanduration,dur - t0)
        t1 = t0 + T

        while (T < 2.0*discrtimestep):
            t0 = random.rand()*dur
            if meanduration == 0:
                meanduration = dur - t0
                
            T = random.rand()*min(meanduration, dur - t0)
            t1 = t0 + T

            if t1 > upperlimit:
                t1 = upperlimit
                if (t1 < t0):
                    temp = t0
                    t0 = t1
                    t1 = temp
                    T = t1 - t0

        # print "\n\nShortcutting iteration", it + 1
        # print t0, t1, t1- t0       
        # interpolate from t0 to t1
        R_beg = lietraj.EvalRotation(t0)
        R_end = lietraj.EvalRotation(t1)
        omega0 = lietraj.EvalOmega(t0)
        omega1 = lietraj.EvalOmega(t1)

        shortcuttraj = lie.InterpolateSO3(R_beg,R_end,omega0,omega1, T)
        #check feasibility only for the new portion

        isincollision = CheckCollisionTraj(robot, shortcuttraj, R_beg, discrtimestep)
        if (not isincollision):
            # a,b,c = lie.ComputeSO3Constraints(shortcuttraj, taumax, discrtimestep)
            abc = TOPPbindings.RunComputeSO3Constraints(str(shortcuttraj),constraintsstring)# discrtimestep)
            a,b,c = lie.Extractabc(abc)

            topp_inst = TOPP.QuadraticConstraints(shortcuttraj, discrtimestep, vmax, list(a), list(b), list(c))
            x = topp_inst.solver
            ret = x.RunComputeProfiles(1,1) 
            if (ret == 1):
                x.resduration
                ## check whether the new one has shorter duration
                if (x.resduration + 0.01 < T): #skip if not shorter than 0.3 s
                    
                    x.ReparameterizeTrajectory()
                    x.WriteResultTrajectory()
                    TOPPed_shortcuttraj = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)

                    newlietraj = ReplaceTrajectorySegment(lietraj,TOPPed_shortcuttraj, t0, t1)  
                    lietraj = newlietraj
                    dur = lietraj.duration
                    #print "*******************************************"
                    print "Success at iteration",it + 1,":", t0, t1,"Deta_t:", t1 - t0 - x.resduration
                    attempt += 1
                    #print "T:", nnotretimable, "; S:", nnotshorter , "; C:", ncollision , "; OK:", attempt
                    #print "*******************************************"
                else:
                    # print "Not shorter"
                    nnotshorter += 1
            else: 
                # print "Not retimable"
                nnotretimable += 1
        else:
            # print "Collision"
            ncollision += 1

            # print "T:", nnotretimable, "; S:", nnotshorter , "; C:", ncollision , "; OK:", attempt
        
    print "\033[1;32mT:", nnotretimable, "; S:", nnotshorter , "; C:", ncollision , "; OK:", attempt, "\033[0m"
    print "\033[1;32m", originalduration - lietraj.duration ,"sec. shorter\033[0m"
    t_sc_end = time.time()
    print "\033[1;32mRunning time:",t_sc_end-t_sc_start, "sec.\033[0m"

    return lietraj

################## REPLACE TRAJECTORY SEGMENT SO3 #############################
def ReplaceTrajectorySegment(originallietraj, trajsegment, t0, t1):
    """ReplaceTrajectorySegment replaces the segment (t0, t1), it returns a LieTraj variable (Rotationlist and Trajectorylist)   """
    assert(t1 > t0)
    newtrajlist = []
    newRlist = []
    i0, rem0 = originallietraj.FindTrajIndex(t0)
    i1, rem1 = originallietraj.FindTrajIndex(t1)

    # print "t" , originallietraj.FindTrajIndex(t0) ##
    # print "t" , originallietraj.FindTrajIndex(t1) ##

    ## check if t0 falls in the first traj. 
    ## if not, insert traj 0 to traj i0 - 1 into newtrajlist
    if i0 > 0:
        for i in range(0,i0):
            newtrajlist.append(originallietraj.trajlist[i])
            newRlist.append(originallietraj.Rlist[i])
    ## remaindertraj0
    newchunkslist = []
    ic0, remc0 = originallietraj.trajlist[i0].FindChunkIndex(rem0)
    # print "c0", originallietraj.trajlist[i0].FindChunkIndex(rem0) ##
          # check if rem0 falls in the first chunk, if not, ...
    if ic0 > 0:
        for c in originallietraj.trajlist[i0].chunkslist[0: ic0]:
            newchunkslist.append(c)
          # remainderchunk0
    remchunk0 = Trajectory.Chunk(remc0, originallietraj.trajlist[i0].chunkslist[ic0].polynomialsvector)
    newchunkslist.append(remchunk0)

    remtraj0 = Trajectory.PiecewisePolynomialTrajectory(newchunkslist) 
    newtrajlist.append(remtraj0)
    newRlist.append(originallietraj.Rlist[i0])


    ## insert trajsegment
    newtrajlist.append(trajsegment)
    newRlist.append(originallietraj.EvalRotation(t0))


######################################
    ## For the traj right after the trajsegment 
    ## remaindertraj1
    newchunkslist = []
    ic1, remc1 = originallietraj.trajlist[i1].FindChunkIndex(rem1)
    newpoly_list = []
    for p in originallietraj.trajlist[i1].chunkslist[ic1].polynomialsvector:
        ## perform variable changing of p(x) = a_n(x)^n + a_(n-1)(x)^(n-1) + ...
        ## by x = y + remc1
        
        a = p.q ## coefficient vector with python convention (highest degree first)
        ## a is a poly1d object
        r = a.r ## polynomial roots
        for i in range(len(r)):
            r[i] = r[i] - remc1
        b = np.poly1d(r, True) ## reconstruct a new polynomial from roots
        ## b is a poly1d object
        b = b*a.coeffs[0] ## multiply back by a_n *** this multiplication does not commute
        
        newpoly = Trajectory.Polynomial(b.coeffs.tolist()[::-1]) ## TOPP convention is weak-term-first
        newpoly_list.append(newpoly)
    remchunk1 = Trajectory.Chunk(originallietraj.trajlist[i1].chunkslist[ic1].duration - remc1, newpoly_list)
    newchunkslist.append(remchunk1)

    ## insert remaining chunk 
    if ic1 < len(originallietraj.trajlist[i1].chunkslist) - 1:
        for c in originallietraj.trajlist[i1].chunkslist[ic1 + 1: len(originallietraj.trajlist[i1].chunkslist)]:
            newchunkslist.append(c)
    ## insert 
    remtraj1 = Trajectory.PiecewisePolynomialTrajectory(newchunkslist)
    newtrajlist.append(remtraj1)
    newRlist.append(originallietraj.Rlist[i1])##ROTATION Should be at originallietraj.Rlist[i1] ##

###############################
    # insert the remainder trajectoris
    if i1 < len(originallietraj.trajlist)-1:
        Rindex = i1+1
        for t in originallietraj.trajlist[i1+1: len(originallietraj.trajlist)]:
            newtrajlist.append(t)
            newRlist.append(originallietraj.Rlist[Rindex])
            Rindex += 1
 
    return lie.LieTraj(newRlist, newtrajlist)





########################### FROM TRAJ LIST TO TRAJSTRING #############################
def TrajStringFromTrajList(trajlist):
    trajectorystring = ""
    for i in range(len(trajlist)):
        trajectorystring += "\n"
        trajectorystring += str(trajlist[i])
    trajectorystring = string.lstrip(trajectorystring) # remove leading "\n"
    return trajectorystring

########################## SAVE LIETRAJ as TEXT_FILEs ################################
## return 2 files: Rlistfilename.txt and trajlistfilename.txt
def SaveLietrajAsTextFiles(lietraj, RlistFilename, trajlistFilename):
    ## Save Rlist
    txtRlist = ""
    for i in range(len(lietraj.Rlist)):
        temp = lietraj.Rlist[i]
        for row in range(0,3):
            separator = ""
            for col in range(0,3):
                txtRlist += separator
                txtRlist += str(temp[row,col])
                separator = " "
            txtRlist += "\n"
    with open(RlistFilename,"wt") as file:
        file.write(txtRlist)
    ## Save trajlist
    txttrajlist = ""
    for i in range(len(lietraj.trajlist)):
        txttrajlist += "t\n"
        txttrajlist += str(lietraj.trajlist[i]) 
        txttrajlist += "\n"
    with open(trajlistFilename,"wt") as file:
        file.write(txttrajlist)
    ## if saved successfully, return true. If not, return false and notify
    return True

######### READ Rlistfilename.txt and trajlistfilename.txt and RETURN a LIETRAJ ######
def ReadLieTrajFiles(Rlistfilename, trajlistfilename):
    ## Read Rlist
    with open(Rlistfilename, 'r') as file:
        data_Rliststring = file.read()
    list = [float(x) for x in data_Rliststring.split()]
    n = len(list)/9
    if (n==0):
        print "\033[91mWRONG DATA, RECHECK ", Rlistfilename, " file!!!\033[0m" 

    Rlist = []
    for i in range(n):
        temp = array([[list[i*9+0], list[i*9+1], list[i*9+2]],
                      [list[i*9+3], list[i*9+4], list[i*9+5]],
                      [list[i*9+6], list[i*9+7], list[i*9+8]]])
        Rlist.append(temp)

    ## Read trajlist
    with open(trajlistfilename, 'r') as file:
        data_trajliststring = file.read()
    buff = StringIO.StringIO(data_trajliststring)
    trajlist = []
    trajstringlist = []
    temptrajstring = ""
    temp = buff.readline() # ingore the first line which contains "t"
    while buff.pos < buff.len:
        temp = buff.readline()
        if (temp != 't\n'):
            temptrajstring += temp
        else:
            trajstringlist.append(temptrajstring)
            temptrajstring = ""
    trajstringlist.append(temptrajstring) # add the last trajstring
    for t_str in trajstringlist:
        traj = Trajectory.PiecewisePolynomialTrajectory.FromString(t_str)
        trajlist.append(traj)

    return  lie.LieTraj(Rlist, trajlist)

################### SAVE SE3 traj########################
## return 2 files: rlistFilename.txt and se3listFilename.txt
def SaveSE3trajAsTextFiles(se3traj, rlist, rlistFilename, se3listFilename):
    ## Save Rlist
    txtrlist = ""
    for i in range(len(rlist)):
        temp = rlist[i]
        for row in range(0,3):
            separator = ""
            for col in range(0,3):
                txtrlist += separator
                txtrlist += str(temp[row,col])
                separator = " "
            txtrlist += "\n"
    with open(rlistFilename,"wt") as file:
        file.write(txtrlist)
    ## Save se3list
    txtse3trajlist = ""
    for i in range(len(se3traj.chunkslist)):
        txtse3trajlist += "t\n"
        txtse3trajlist += str(se3traj.chunkslist[i]) 
        txtse3trajlist += "\n"
    with open(se3listFilename,"wt") as file:
        file.write(txtse3trajlist)
    ## if saved successfully, return true. If not, return false and notify
    return True


######### READ rlistFilename.txt and se3trajFilename.txt and RETURN a SE3TRAJ ######
def ReadSE3TrajFiles(rlistfilename, se3trajfilename):
    ## Read rlist
    with open(rlistfilename, 'r') as file:
        data_rliststring = file.read()
    list = [float(x) for x in data_rliststring.split()]
    n = len(list)/9
    if (n==0):
        print "\033[91mWRONG DATA, RECHECK ", rlistfilename, " file!!!\033[0m" 

    rlist = []
    for i in range(n):
        temp = array([[list[i*9+0], list[i*9+1], list[i*9+2]],
                      [list[i*9+3], list[i*9+4], list[i*9+5]],
                      [list[i*9+6], list[i*9+7], list[i*9+8]]])
        rlist.append(temp)

    ## Read trajlist
    with open(se3trajfilename, 'r') as file:
        data_se3trajliststring = file.read()
    buff = StringIO.StringIO(data_se3trajliststring)
    se3trajlist = []
    se3trajstringlist = []
    tempse3trajstring = ""
    temp = buff.readline() # ingore the first line which contains "t"
    while buff.pos < buff.len:
        temp = buff.readline()
        if (temp != 't\n'):
            tempse3trajstring += temp
        else:
            se3trajstringlist.append(tempse3trajstring)
            tempse3trajstring = ""
    se3trajstringlist.append(tempse3trajstring) # add the last trajstring
    for t_str in se3trajstringlist:
        se3traj = Trajectory.PiecewisePolynomialTrajectory.FromString(t_str)
        se3trajlist.append(se3traj)
    se3traj = Trajectory.PiecewisePolynomialTrajectory.FromString(TrajStringFromTrajList(se3trajlist))
    return  se3traj, rlist

#########################PLOT SE3 ###################################
def PlotSE3(se3traj, rlist,  dt = 0.01, figstart=0,vmax=[],accelmax=[],taumax=[],fmax=[], inertia = None, m = None):
    transtraj, rottraj = TransRotTrajFromSE3Traj(se3traj)
    lietraj = lie.SplitTraj2(rlist, rottraj)
    
    lietraj.Plot(dt,figstart,vmax[:3],accelmax,taumax,inertia)
    
    figure(figstart+3)
    clf()
    tvect = arange(0, transtraj.duration + dt, dt)
    qdvect = array([transtraj.Evald(t) for t in tvect])
    plt.plot(tvect, qdvect[:,0], '--', label = r'$v^1$',linewidth=2)
    plt.plot(tvect, qdvect[:,1], '-.', label = r'$v^2$',linewidth=2)
    plt.plot(tvect, qdvect[:,2], '-', label = r'$v^3$',linewidth=2)
    plt.legend()
    ylabel('Translation velocities (m/s)')
    xlabel('Time (s)')
    for v in vmax[:3]:
        plt.plot([0, transtraj.duration],[v, v], '-.',color = 'k')
    for v in vmax[:3]:
        plt.plot([0, transtraj.duration],[-v, -v], '-.',color = 'k')

    figure(figstart+4)
    clf()
    qddvect = array([transtraj.Evaldd(t) for t in tvect])
    plt.plot(tvect, qddvect[:,0], '--', label = r'$f^1$',linewidth=2)
    plt.plot(tvect, qddvect[:,1], '-.', label = r'$f^2$',linewidth=2)
    plt.plot(tvect, qddvect[:,2], '-', label = r'$f^3$',linewidth=2)
    plt.legend()
    ylabel('Forces (N)')
    xlabel('Time (s)')
    for v in fmax[:3]:
        plt.plot([0, transtraj.duration],[v, v], '-.',color = 'k')
    for v in fmax[:3]:
        plt.plot([0, transtraj.duration],[-v, -v], '-.',color = 'k')

###########################CheckIntersection ##########################
def CheckIntersection(interval0, interval1):
    """CheckIntersection checks whether interval0 intersects interval1.
    """
    
    if (np.max(interval0) < np.min(interval1)):
        return False

    elif (np.max(interval1) < np.min(interval0)):
        return False
    
    else:
        return True
