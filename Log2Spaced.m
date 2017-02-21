function X = Log2Spaced(Start, Stop, N)
%Returns logarithmically-spaced points from 2^Start to 2^Stop.
%
%inputs:
%Start - (scalar) Starting power of two.
%Stop - (scalar) Stopping power of two.
%N - (scalar) Number of points.
%
%outputs:
%X - (N-length float) log-spaced points between 2^Start and 2^Stop.
%
%Author: Lee Cooper, Emory University.


X = (2).^ [Start+(0:N-2)*(Stop-Start)/(N-1), Stop];
