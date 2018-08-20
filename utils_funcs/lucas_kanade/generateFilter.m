function filter = generateFilter (alpha)
% Name: 		generateFilter
% Purpose: 	Get 1D filtering for pyramid generation
% Returns:  Desired 1D filter
% Parameters: 
%           alpha: the only parameter determining the filter shape
% Notes: 
% 				Here we use only 1D filtering, because the 2D filter is separable
% 				By doing 2 1D filtering we can significantly reduce the time complexity
% 				without portability and understandability loss.

filter = [0.25-alpha/2 0.25 alpha 0.25 0.25-alpha/2];