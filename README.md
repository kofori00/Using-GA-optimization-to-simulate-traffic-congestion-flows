# Simulating-traffic-congestion-flows
Traffic congestion on freeways is a critical problem due to its negative impact on the
environment and many other important consequences (higher delays, waste of fuel, higher
accident risk probability, etc.). Since the construction of new freeways is not always a viable
option or is too costly, other solutions have to be found. In many cases, the use of dynamic
trac control measures such as ramp metering, variable speed limits, reversible lanes, and
route guidance may be a cost-effcient and eective solution.
In general, dynamic traffic control uses measurements of the traffic conditions over time
and computes dynamic control signals to in
uence the behavior of the drivers and to generate
a response in such a way that the performance of the network is improved, by reducing delays,
emissions, fuel consumption etc.

## Description 
The most used control input in freeway traffic control is ramp metering, a device (usually
a basic traffic light) located at an on-ramp that regulates the 
ow of traffic entering the
freeway according to the current traffic conditions. Ramp metering systems have proved to
be successful in decreasing trac congestion and improving driver safety.
Another promising control measure in freeway traffic control is the use of Variable Speed Limits (VSLs), overhead signs showing a varying speed limit that changes according to the current road conditions.

## Goal 
The goal of this assignment is to compute the values of the speed limits and the ramp
metering rates (i.e. the percentage of ramp 
ow that is allowed to enter) that minimize the
congestion that appears on a freeway stretch.

## Model 
The macroscopic model METANET will be used to model the behavior of the freeway.
METANET represents the traffic network as a graph where the segments are indexed by an
index i and have a length of L meters with $\lamba$ lanes.   

