function [ J, theta, Jx, Jy ] = findDerivatives( I, Gx, Gy )
%% findDerivative 
%  Calculate the x and y derivatives of image I
%  Input: 
%   I - a gray-scale image; 
%   Gx, Gy - Gaussian filter in x and y direction
%  Output:
%   J - magnitute of derivative
%   theta - orientation of derivative
%   Jx, Jy - magnitute of derivative along x and y axis
%
%  Jiani Li, 09/27/2016

% Normalize to Gaussian kernels
Gx = Gx/sum(Gx); Gy = Gy/sum(Gy);

dGx = gradient(Gx);
dGy = gradient(Gy);

Jx = imfilter(double(I),Gy,'conv','replicate');
Jx = imfilter(Jx,dGx,'conv','replicate');
Jy = imfilter(double(I),Gx,'conv','replicate');
Jy = imfilter(Jy,dGy,'conv','replicate');


J = sqrt(Jx.*Jx + Jy.*Jy);
theta = atan2(Jy,Jx);


end

