function [blk,A,C,b] = sdpa_to_vsdp(fname)
% SDPA_TO_VSDP  Import problem in SDPA sparse format to VSDP's format.
%
%   [blk,A,C,b] = sdpa_to_vsdp(fname)  The input 'fname' name of the file
%      containing SDP data in SDPA format.  The problem data (blk,A,C,b) of
%      VSDP's block-diagonal problem is described in 'mysdps.m'.
%
%   Example:
%
%       fname = fullfile('sdplib_vsdp', 'arch4.dat-s');
%       [blk,A,C,b] = sdpa_to_vsdp(fname);
%
%   See also mysdps.

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)

%Open the file for input
compressed = 0;
if exist(fname,'file')
  fid = fopen(fname,'r');
elseif exist([fname,'.Z'],'file')
  compressed = 1;
  unix(['uncompress ',fname,'.Z']);
  fid = fopen(fname,'r');
elseif exist([fname,'.gz'],'file')
  compressed = 2;
  unix(['gunzip ',fname,'.gz']);
  fid = fopen(fname,'r');
else
  fprintf('** Problem not found, please specify the correct path. \n');
  return;
end

%Clean up special characters and comments from the file

datavec = fscanf(fid,'%c');
linefeeds = findstr(datavec,char(10));
comment_chars = '*"=';
cumidx = [];
for i=1:length(comment_chars)
  idx = findstr(datavec,comment_chars(i));
  cumidx = [cumidx,idx];
end
for j=length(cumidx):-1:1
  if (cumidx(j)==1) || (strcmp(datavec(cumidx(j)-1),char(10)))
    datavec(cumidx(j):linefeeds(min(find(cumidx(j)<linefeeds))))='';
  else
    datavec(cumidx(j):linefeeds(min(find(cumidx(j)<linefeeds)))-1)='';
  end
end
special_chars = ',{}()';
cumidx=[];
for i=1:length(special_chars)
  idx = findstr(datavec,special_chars(i));
  cumidx = [cumidx,idx];
end
datavec(cumidx) = blanks(length(cumidx));
clear linefeeds;

%Close the file

fclose('all');
if compressed==1; unix(['compress ',fname]); end
if compressed==2; unix(['gzip ',fname]); end

%Next, read in basic problem size parameters.

datavec = sscanf(datavec,'%f');
if size(datavec,1) < size(datavec,2); datavec = datavec'; end
m = datavec(1);                        %the number of dual variables
numblk  = datavec(2);                  %the number of blocks in a matrix
blksize = abs(datavec(2+(1:numblk)));       %the block structure vector
if size(blksize,1) > size(blksize,2); blksize = blksize'; end
A=cell(numblk,m);
C=cell(numblk,1);

%Get input  b.
idxstrt = 2+numblk;
b = datavec(idxstrt+(1:m));
idxstrt = idxstrt+m;
b = -b;

%Construct blk
%
denumblk = length(blksize);
for p = 1:denumblk
  sj = blksize(p);
  blk{p,1} = 's'; blk{p,2} = sj;
end
%blk{1,1} = 's'; blk{1,2} = blksize;
%Construct A and C
len = length(datavec);
Y = reshape(datavec(idxstrt+1:len),5,(len-idxstrt)/5)';
clear datavec;
Y = sortrows(Y,[1 2]);
matidx = [0; find(diff(Y(:,1)) ~= 0); size(Y,1)];

for k = 1:length(matidx)-1
  idx = (matidx(k)+1 : matidx(k+1));
  Ytmp  = Y(idx,1:5);
  matno = Ytmp(1,1);
  Ytmp2 = Ytmp(:,2);
  for p = 1:numblk
    n  = blksize(p);
    idx = find(Ytmp2 == p);
    ii = Ytmp(idx,3); jj = Ytmp(idx,4); vv =Ytmp(idx,5);
    
    idxtmp = find(ii > jj);
    if ~isempty(idxtmp)
      tmp = jj(idxtmp);
      jj(idxtmp) = ii(idxtmp); ii(idxtmp) = tmp;
    end
    tmp = -sparse(ii,jj,vv,n,n);
    tmp = tmp + triu(tmp,1)';
    
    if (matno == 0)
      C{p,1} = tmp;
    else
      % A{p,1}(:,matno) = tmp;
      A{p,matno} = tmp;
    end
  end
end

A = A';

end
