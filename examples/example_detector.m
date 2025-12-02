function [idx] = QRSDetect(fileName,m, normCnst)
  S = load(fileName);
  sigO = S.val(1,:).* 5;
  sigLen = size(sigO,2);
  sig1 = zeros (1,sigLen);
  sig2 = zeros (1,sigLen);
  sig3 = zeros (1,sigLen);
  sig4 = zeros (1,sigLen);
  dmax=0;
  idx=[];
  %1
  for i=2:sigLen-1
      sig1(1,i)=(sigO(1,i-1)+2*sigO(1,i)+sigO(1,i+1))/4.0;
  end
  %2
  for i=m+1:sigLen-m
    tmpS=0;
    for k=i-m:i+m tmpS=tmpS+sig1(1,k); end
    tmpS=tmpS/(2*m+1);
    sig2(1,i)=(sig1(1,i)-tmpS)*(sig1(1,i)-tmpS);
  end
  %3
  for i=m+1:sigLen-m
    tmpS=0;
    for k=i-m:i+m tmpS = tmpS + sig2(1,k); end
    tmpS=tmpS/64.0;
    sig3(1,i)=sig2(1,i)*tmpS*tmpS/64.0;
  end
  %4
  for i=m+1:sigLen-m
    if ((sig1(1,i)-sig1(1,i-m))*(sig1(1,i)-sig1(1,i+m))>=0) sig4(1,i)=sig3(1,i); end
    if (sig4(1,i)>dmax) dmax=sig4(1,i); end
  end
  dmax=dmax/normCnst;
  for i=1:sigLen if sig4(1,i)>dmax sig4(1,i)=1; else sig4(1,i)=0; end 
  end

  qrsInt=36;
  i=1;
  while i<=sigLen
    if (sig4(1,i)>0)
      for j=i:i+qrsInt-1 if (j<=sigLen) sig4(1,j)=1.0; end 
      end
      max=i+1;
      for k=i+1:i+qrsInt-1 
          if (k<=sigLen) 
              if ((sigO(1,max)-sigO(1,max-1))<(sigO(1,k)-sigO(1,k-1))) max=k;  end 
          end
      end
      for k=i:i+qrsInt-1 if (k<=sigLen) sig4(1,k)=0; end; end
      sig4(1,max)=1;
      i=i+qrsInt;
    else
        i=i+1;
    end
  end %while i=1
  for i=1:sigLen if (sig4(1,i)>0) idx=[idx i-1];end; end
end