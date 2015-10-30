# Interpolation in SO(3) following Park and Ravani
import time

import TOPP
from TOPP import Trajectory
import bisect
from pylab import *
from numpy import *
import pdb
import matplotlib.pyplot as plt




Eps = zeros((3,3,3))
Eps[0,2,1] = -1
Eps[0,1,2] = 1
Eps[1,0,2] = -1
Eps[1,2,0] = 1
Eps[2,1,0] = -1
Eps[2,0,1] = 1


class LieTraj():
    def __init__(self, Rlist,trajlist):
        self.Rlist = Rlist
        self.trajlist = trajlist
        self.duration = 0
        self.trajcumulateddurationslist = []
        for t in trajlist:
            self.trajcumulateddurationslist.append(self.duration)
            self.duration += t.duration

    def FindTrajIndex(self, s):
        if s == 0:
            s = 1e-10
        i = bisect.bisect_left(self.trajcumulateddurationslist, s) - 1
        remainder = s - self.trajcumulateddurationslist[i]
        return i, remainder

    # Rotation
    def EvalRotation(self,s):
        i, remainder = self.FindTrajIndex(s)
        return(dot(self.Rlist[i],expmat(self.trajlist[i].Eval(remainder))))

    # Velocity in body frame
    def EvalOmega(self,s):
        i, remainder = self.FindTrajIndex(s)
        r = self.trajlist[i].Eval(remainder)
        rd = self.trajlist[i].Evald(remainder)
        return dot(Amat(r),rd)

    # Acceleration in body frame
    def EvalAlpha(self,s):
        i, remainder = self.FindTrajIndex(s)
        r = self.trajlist[i].Eval(remainder)
        rd = self.trajlist[i].Evald(remainder)
        rdd = self.trajlist[i].Evaldd(remainder)
        return dot(Bmat(r),rdd) + dot(rd,tensordot(Ctensor(r),rd,([2],[0])))

    # Torques
    def EvalTorques(self,s,I):
        i, remainder = self.FindTrajIndex(s)
        r = self.trajlist[i].Eval(remainder)
        rd = self.trajlist[i].Evald(remainder)
        rdd = self.trajlist[i].Evaldd(remainder)
        omega = dot(Amat(r),rd)
        alpha =  dot(Bmat(r),rdd) + dot(rd,tensordot(Ctensor(r),rd,([2],[0])))
        return dot(I,alpha) + cross(omega,dot(I,omega))

    
    def Plot(self,dt=0.01,figstart=0,vmax=[],accelmax=[],taumax=[],I=None):

        tvect = arange(0, self.duration + dt, dt)
        omegavect = array([self.EvalOmega(t) for t in tvect])
        figure(figstart)
        clf()
        
        plt.plot(tvect,omegavect[:,0],'--',label = '$\omega^1$',linewidth = 2)
        plt.plot(tvect,omegavect[:,1],'-.',label = '$\omega^2$',linewidth = 2)
        plt.plot(tvect,omegavect[:,2],'-',label ='$\omega^3$',linewidth = 2)
        plt.legend()
        for v in vmax:
            plt.plot([0, self.duration],[v, v], '-.',color = 'k')
        for v in vmax:
            plt.plot([0, self.duration],[-v, -v], '-.',color = 'k')
        ylabel('Angular velocities (rad/s)')
        xlabel('Time (s)')

        alphavect = array([self.EvalAlpha(t) for t in tvect])
        figure(figstart+1)
        clf()
        plt.plot(tvect,alphavect[:,0],'--',label = '$\dot \omega^1$',linewidth = 2)
        plt.plot(tvect,alphavect[:,1],'-.',label = '$\dot \omega^2$',linewidth = 2)
        plt.plot(tvect,alphavect[:,2],'-',label = '$\dot \omega^3$',linewidth = 2)      
        plt.legend()
        for a in accelmax:
            plt.plot([0, self.duration],[a, a], '-.',color = 'k')
        for a in accelmax:
            plt.plot([0, self.duration],[-a, -a], '-.',color = 'k')
        ylabel('Angular accelerations (rad/s^2)')
        xlabel('Time (s)')

        if not(I is None):
            torquesvect = array([self.EvalTorques(t,I) for t in tvect])
            figure(figstart+2)
            clf()
            plt.plot(tvect,torquesvect[:,0],'--',label = r'$\tau^1$',linewidth = 2)
            plt.plot(tvect,torquesvect[:,1],'-.',label = r'$\tau^2$',linewidth = 2)
            plt.plot(tvect,torquesvect[:,2],'-',label = r'$\tau^3$',linewidth = 2)
            plt.legend()
            
            for tau in taumax:
                plt.plot([0, self.duration],[tau, tau], '-.',color = 'k')
            for tau in taumax:
                plt.plot([0, self.duration],[-tau, -tau], '-.',color = 'k')
            ylabel('Torques (N.m)')
            xlabel('Time (s)')
        
def SplitTraj(Rlist,traj):
    trajlist = []
    chunkindex = 0
    clist = []
    for i in range(len(Rlist)-1):
        while chunkindex <  len(traj.chunkslist):
            chunkcur = traj.chunkslist[chunkindex]
            chunknext = traj.chunkslist[chunkindex+1]
            clist.append(chunkcur)
            chunkindex += 1
            if(norm(dot(Rlist[i],expmat(chunkcur.Eval(chunkcur.duration)))-dot(Rlist[i+1],expmat(chunknext.Eval(0)))))< 1e-8:
                trajlist.append(Trajectory.PiecewisePolynomialTrajectory(clist))
                clist = []
                break
    # Last traj
    clist = []
    while chunkindex < len(traj.chunkslist):
        clist.append(traj.chunkslist[chunkindex])
        chunkindex += 1
    trajlist.append(Trajectory.PiecewisePolynomialTrajectory(clist))
        
    return LieTraj(Rlist,trajlist)
      
def SplitTraj2(Rlist,traj): 
    trajlist = []
    chunkindex = 0
    clist = []
    for i in range(len(Rlist)-1):
        while chunkindex <  len(traj.chunkslist):
            chunkcur = traj.chunkslist[chunkindex]
            chunknext = traj.chunkslist[chunkindex+1]
            
            clist.append(chunkcur)
            chunkindex += 1
            if(norm(dot(Rlist[i],expmat(chunkcur.Eval(chunkcur.duration)))-dot(Rlist[i+1],expmat(chunknext.Eval(0)))))< 1e-1:
                trajlist.append(Trajectory.PiecewisePolynomialTrajectory(clist))
                clist = []
                break
    # Last traj
    clist = []
    while chunkindex < len(traj.chunkslist):
        clist.append(traj.chunkslist[chunkindex])
        chunkindex += 1
    trajlist.append(Trajectory.PiecewisePolynomialTrajectory(clist))
        
    return LieTraj(Rlist,trajlist)


def skewfromvect(r):
    return array([[0,-r[2],r[1]],[r[2],0,-r[0]],[-r[1],r[0],0]])

def vectfromskew(R):
    return array([R[2,1],R[0,2],R[1,0]])

def expmat(r):
    nr = linalg.norm(r)
    if(nr<=1e-10):
        return eye(3)
    R = skewfromvect(r)
    return eye(3) + sin(nr)/nr*R + (1-cos(nr))/(nr*nr)*dot(R,R)

def logvect(R):
    if(abs(trace(R)+1)>1e-10):
        if(linalg.norm(R-eye(3))<=1e-10):
            return zeros(3)
        else:
            phi = arccos((trace(R)-1)/2)
            return vectfromskew(phi/(2*sin(phi))*(R-R.T))
    else:
        eigval, eigvect = linalg.eig(R)
        for (i,val) in enumerate(eigval):
            if abs((val-1)) <= 1e-10:
                return pi*real(eigvect[:,i])
                
def Amat(r):
    nr = linalg.norm(r)
    if(nr<=1e-10):
        return eye(3)
    R = skewfromvect(r)
    return eye(3) - (1-cos(nr))/(nr*nr)*R + (nr-sin(nr))/(nr*nr*nr)*dot(R,R)

def Bmat0(r):
    nr = linalg.norm(r)
    R = skewfromvect(r)
    return eye(3) + (1-cos(nr))/(nr*nr)*R + (nr-sin(nr))/(nr*nr*nr)*dot(R,R)

def Bmat(r):
    nr = linalg.norm(r)
    R = skewfromvect(r)
    return eye(3) - (1-cos(nr))/(nr*nr)*R + (nr-sin(nr))/(nr*nr*nr)*dot(R,R)

def Ctensor(r):
    nr = linalg.norm(r)
    nr2 = nr*nr
    nr3 = nr2*nr
    nr4 = nr3*nr
    nr5 = nr4*nr
    R = skewfromvect(r)
    C1 = -(nr-sin(nr))/nr3 * dot(Eps,R)
    C2 = -(2*cos(nr)+nr*sin(nr)-2)/nr4 * TensorProd(r,R)
    C3 = (3*sin(nr)-nr*cos(nr) - 2*nr)/nr5 * TensorProd(r,dot(R,R))
    return C1+C2+C3

def Cterm(r,rd):
    nr = linalg.norm(r)
    nr2 = nr*nr
    nr3 = nr2*nr
    nr4 = nr3*nr
    nr5 = nr4*nr
    C1 = (nr-sin(nr))/nr3 * cross(rd,cross(r,rd))
    C2 = -(2*cos(nr)+nr*sin(nr)-2)/nr4 * dot(r,rd)*cross(r,rd)
    C3 = (3*sin(nr)-nr*cos(nr) - 2*nr)/nr5 * dot(r,rd)*cross(r,cross(r,rd))
    return C1+C2+C3
   

def omega(r,rd):
    return dot(Amat(r),rd)

def alpha(r,rd,rdd):
    return dot(Bmat(r),rdd) + dot(rd,tensordot(Ctensor(r),rd,([2],[0])))

def tau(r,rd,rdd,I):
    omega0 = omega(r,rd)
    return dot(I,alpha(r,rd,rdd)) + cross(omega0,dot(I,omega0))


def InterpolateSO3(R0,R1,omega0,omega1,T):

    r1 = logvect(dot(R0.T,R1))
    u = linalg.solve(Amat(r1),omega1*T)

    c = omega0*T
    M = array([[1,0,0,1,0,0],
               [0,1,0,0,1,0],
               [0,0,1,0,0,1],
               [3,0,0,2,0,0],
               [0,3,0,0,2,0],
               [0,0,3,0,0,2]])
    y = array([r1[0]-c[0],
               r1[1]-c[1],
               r1[2]-c[2],
               u[0]-c[0],
               u[1]-c[1],
               u[2]-c[2]])

    x = linalg.solve(M,y)
    a = x[:3]
    b = x[3:]
    T2 = T*T
    T3 = T2*T
    polylist = []

    for i in range(3):
        polylist.append(Trajectory.Polynomial([0,c[i]/T,b[i]/T2,a[i]/T3]))
    chunk = Trajectory.Chunk(T,polylist)
    return Trajectory.PiecewisePolynomialTrajectory([chunk])
    


def EvalRotation(R0,traj,t):
    return(dot(R0,expmat(traj.Eval(t))))
    

def TensorProd(a,A):
    res = zeros((3,3,3))
    for i in range(3):
        res[i,:,:] = a[i]*A
    return res


def ComputeSO3Torques(rtraj, I = None, dt=0.01):
    if I is None:
        I = eye(3)
    tvect = arange(0, rtraj.duration + dt, dt)
    tauvect = []
    for t in tvect:
        r = rtraj.Eval(t)
        rd = rtraj.Evald(t)
        rdd = rtraj.Evaldd(t)
        taut = tau(r,rd,rdd,I)
        tauvect.append(taut)
    return tvect,array(tauvect)
    

def ComputeSO3Constraints(rtraj, taumax, discrtimestep, I = None):
    ndiscrsteps = int((rtraj.duration + 1e-10) / discrtimestep) + 1
    a = zeros((ndiscrsteps,6))
    b = zeros((ndiscrsteps,6))
    c = zeros((ndiscrsteps,6))
    for i in range(ndiscrsteps):
        t = i * discrtimestep
        r = rtraj.Eval(t)
        rd = rtraj.Evald(t)
        rdd = rtraj.Evaldd(t)
        nr = linalg.norm(r)
        nr2 = nr*nr
        nr3 = nr2*nr
        nr4 = nr3*nr
        nr5 = nr4*nr
        R = skewfromvect(r)

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
        
        a[i,:3] = at
        a[i,3:] = -at
        b[i,:3] = bt
        b[i,3:] = -bt
        c[i,:3] = -taumax
        c[i,3:] = -taumax
    return a, b, c

def RandomQuat():
    s = random.rand()
    sigma1 = sqrt(1-s)
    sigma2 = sqrt(s)
    theta1 = 2*pi*random.rand()
    theta2 = 2*pi*random.rand()
    w = cos(theta2)*sigma2
    x = sin(theta1)*sigma1
    y = cos(theta1)*sigma1
    z = sin(theta2)*sigma2
    return array([w,x,y,z])
    

def InterpolateSO3ZeroOmega(R0,R1,T):
    r = logvect(dot(R0.T,R1))
    a = ones(3)*(-2)
    b = ones(3)*3
    T2 = T*T
    T3 = T2*T
    polylist = []
    for i in range(3):
        polylist.append(Trajectory.Polynomial([0,0,r[i]*b[i]/T2,r[i]*a[i]/T3]))
    chunk = Trajectory.Chunk(T,polylist)
    return Trajectory.PiecewisePolynomialTrajectory([chunk])

def Extractabc(abc):
    lista = [float(x) for x in abc[0].split()]
    listb = [float(x) for x in abc[1].split()]
    listc = [float(x) for x in abc[2].split()]
    n= len(lista)/6
    a = zeros((n,6))
    b = zeros((n,6))
    c = zeros((n,6))
    for i in range(n):
        a[i,:] = lista[i*6:i*6+6]
        b[i,:] = listb[i*6:i*6+6]
        c[i,:] = listc[i*6:i*6+6]
    return a, b, c
