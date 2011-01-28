function [ pdf_total ] = nfoldconvolution( pdf_matrix )
%nfoldconvolution Perform a convolution sum of PDF
%   Given a matrix of probability distribution functions, this function
%   performs a convolution sum and outputs a single combined probability
%   distribution function.

s = size(pdf_matrix,1);

pdf_total = conv(pdf_matrix(1,:),pdf_matrix(2,:));

for i=3:s,
    pdf_total = conv(pdf_total,pdf_matrix(i,:));
end

end

