%% this function is to produce the full array
function Fa=FullAarry(basicE,k)
% basicE is the number of basic element
% k is the length of the produced string
%n=length(basicE);
substring=basicE;
for i=1:k-1
    substring=SubString(substring,basicE);
end
Fa=substring;
end

function substring=SubString(sub,basicE)
N1=length(sub);
N2=length(basicE);
t=1;
for i=1:N1
    for j=1:N2
        substring{t,1}=strcat(sub{i,1},basicE{j,1});
        t=t+1;
    end
end
end
    