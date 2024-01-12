

disp(' ')
disp('Simplifying Block Model...');
disp(' ')

segnds=[159 204;
        159 184;
        166 219;
        231 244;
        223 244;
        91 109;
        99 115;
        4  28
        71 119
%         51 98
        66 71
        54 58
        28 61
        28 79
        255 266
        170 120
        138 170
        144 162
        119 106
        202 209
        231 238
        174 224
        178 183
        124 159
        217 220
        18 31
        92 104
        47 76
        11 41
        48 50
        84 89
        93 131
        159 201
        210 213
        131 171
        171 176
        180 203
        131 177
        201 204];

for i=1:size(segnds,1)
    [blocksnew,faultsnew,blocksout,faultsout]=RemoveBlock(blocks,nodes,faults,segnds(i,:)');
    blocks=blocksnew;
    faults=faultsnew;
end

