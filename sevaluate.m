
function [Jaccard,Dice,rfp,rfn]=sevaluate(m,o)
m=m(:);
o=o(:);
common=sum(m & o); 
union=sum(m | o); 
cm=sum(m); 
co=sum(o); 
Jaccard=common/union;
Dice=(2*common)/(cm+co);
rfp=(co-common)/cm;
rfn=(cm-common)/cm;