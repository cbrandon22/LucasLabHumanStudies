% function to do gaussian smoothing of each column of indata data with 
% given sigma. Periodic boundary condition, end wrap around. Faster version.
% 130398 [outdata] = fgsmooth(indata,sigma) global GAUSS ABINS works.
%
	function [outdata] = fgsmooth(indata,sigma)
%
	[GAUSS] = fgauss(sigma);

  	[len,wid] = size(indata);
	if len == 1
		indata = indata';
  		[len,wid] = size(indata);
	end

	posax = 1:len;
	negax = posax;

	outdata = indata*GAUSS(1);
	for k = 2:length(GAUSS)
		posax = [posax(2:len) posax(1)];
		negax = [negax(len) negax(1:len-1)];
		outdata = outdata + (indata(posax,:) + indata(negax,:))*GAUSS(k);
	end

