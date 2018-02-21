% For given sigma, it generates normalized gaussian distbn of length 3*sigma+1
% [GAUSS]=fgauss(sigma) final tested version on 04Sept96. again modified on
% 130398 when GAUSS and ABINS are made local and normalizn constnat is made
% automatic, thanks to the Herminator's suggestion.

	function [GAUSS] = fgauss(sigma)
%
  ABINS = round(3*sigma+1);
  if sigma <= 0
    GAUSS = ones([1,1:ABINS]);
    return;
  end
  for i = 0:ABINS+1
    GAUSS(i+1) = exp(-(i*i)/(sigma*sigma));
  end
  GAUSS = GAUSS/(2*sum(GAUSS(2:ABINS+1))+GAUSS(1));
