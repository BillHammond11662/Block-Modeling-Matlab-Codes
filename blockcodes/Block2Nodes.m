function [nds]=Block2Nodes(blocks,bnum)
%
% nds = Block2Nodes(blocks,bnum);

blocknames = fieldnames(blocks);
%eval(['nds=blocks.' char(blocknames(bnum)) ';']);
nds = blocks.(char(blocknames(bnum)));