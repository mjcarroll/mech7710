function [ mean, var, c_moment, mean_sq ] = statistics( pdf )
s = length(pdf);
mean = 0;
mean_sq = 0;
var = 0;
c_moment = 0;

for i=1:length(pdf),
    mean = mean + i*pdf(i);
    mean_sq = mean_sq + (i^2) * pdf(i); 
end
for i=1:length(pdf),
    c_moment = c_moment + (i - mean) * pdf(i);
    var = var + (i - mean)^2 * pdf(i);
end
end