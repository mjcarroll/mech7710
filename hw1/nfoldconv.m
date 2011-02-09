function [ pdf_total ] = nfoldconv( pdf_matrix )
%NFOLDCONV Perform a convolution sum of PDF functions.
%   pdf_total = NFOLDCONV(pdf_matrix) iterates over a matrix 
%   of PDFs and convolves each into a single PDF
%   representing the sum of the random variables.

s = size(pdf_matrix,1);

pdf_total = conv(pdf_matrix(1,:),pdf_matrix(2,:));

for i=3:s,
    pdf_total = conv(pdf_total,pdf_matrix(i,:));
end

% Correct MATLAB's convolution function chopping leading zeros off.
pdf_total = [zeros(1,s-1),pdf_total];

assert(sum(pdf_total) - 1.0 <= eps, ...
    'Combined PDF has a sum greater than 1, recheck PDF matrix');

end

