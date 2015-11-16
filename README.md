TOPP FOR RIGID-BODY MOTIONS
------------

This library includes tools to plan fast trajectories on the the space of rotation matrices SO(3) and the space of three-dimensional rigid body motions SE(3) under kinodynamic constraints (bounds on angular velocities, accelerations and torques) in a cluttered environment.

For further details, please see more at our paper <<link to paper>>

Requirements
------------

- As this implementation is an extention of TOPP (time-optimal Path Parameterization), read instructions in following link to install TOPP, OpenRAVE and prerequisites:
https://github.com/quangounet/TOPP

- Clone this TOPP-SO3 folder

$ sudo python setup.py install


Examples
------------
Please try the test files in the examples folder (test-SO3.py and test-SE3.py)

Below is the video demonstrating the resulting trajectories found by our algorithm.

https://youtu.be/heM7uxGrfVc
