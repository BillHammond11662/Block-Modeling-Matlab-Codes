function blocksout = NiceBlockNames(blocks)
% blocksout = NiceBlockNames(blocks)
%
% gives nice short and neat names to blocks

blocksout = struct([]);
bnamesold = fieldnames(blocks);
for i=1:length(bnamesold)
   bnamenew = ['block' sprintf('%.0f',i)];
   ndlist = blocks.(bnamesold{i});
   if isempty(blocksout)
      blocksout = struct(bnamenew,ndlist);
   else
      blocksout.(bnamenew) = ndlist;
   end
end
